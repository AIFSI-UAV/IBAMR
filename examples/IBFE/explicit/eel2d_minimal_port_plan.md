# IBFE 版本 eel2d 复现：最小改造清单（文件级别 + 关键函数骨架）

> 目标：在 **尽量少改动** 的前提下，用 `IBFEMethod + IBFEDirectForcingKinematics` 复现 `examples/ConstraintIB/eel2d` 的游动轨迹与尾迹（先做运动学复现，不追求弹性肌肉机理）。

---

## 0. 复现策略（一句话）

- 保留 `ConstraintIB/eel2d` 的波形表达（`body_shape_equation` 和 `deformation_velocity_function_*`）思想。
- 在 IBFE 中改为：
  1) 用 FE 网格表示鱼体；
  2) 用 `IBFEDirectForcingKinematics` 施加刚体自由度；
  3) 用一个自定义“目标形变/速度”更新节点参考位置（或通过体力项逼近）。

最小可行版本优先采用“**直接给节点目标速度/位置**”路线，避免先引入复杂本构。

---

## 1. 文件级最小改造清单

## 1.1 新增示例目录（建议）

- 新建：`examples/IBFE/explicit/eel2d_dfk/`
  - `example.cpp`
  - `Eel2dKinematics.h`
  - `Eel2dKinematics.cpp`
  - `input2d`
  - `eel2d.exo`（或你现有的二维 FE 网格文件）
  - `CMakeLists.txt`

> 如果你想更“最小”，也可以基于 `examples/IBFE/explicit/ex11` 直接改名复制，保留其工程组织方式。

## 1.2 构建系统

- `examples/IBFE/explicit/CMakeLists.txt`
  - 增加子目录：`add_subdirectory(eel2d_dfk)`。

- `examples/IBFE/explicit/eel2d_dfk/CMakeLists.txt`
  - 参照 `ex11` 链接同类库。

---

## 2. 关键代码骨架

## 2.1 `example.cpp`：主程序骨架

```cpp
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBFEDirectForcingKinematics.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
// ... 其他常规 include
#include "Eel2dKinematics.h"

int main(int argc, char* argv[])
{
    // 1) 初始化 AppInitializer / 网格层级 / NS 积分器

    // 2) 读取 FE mesh（二维）
    auto mesh = ...; // libMesh mesh 初始化

    Pointer<IBFEMethod> ib_method_ops = new IBFEMethod(
        "IBFEMethod",
        app_initializer->getComponentDatabase("IBFEMethod"),
        mesh,
        max_levels);

    Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
        "IBHierarchyIntegrator",
        app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
        ib_method_ops,
        navier_stokes_integrator);

    // 3) Direct Forcing Kinematics
    Pointer<IBFEDirectForcingKinematics> dfk = new IBFEDirectForcingKinematics(
        "eel_dfk",
        app_initializer->getComponentDatabase("EelIBFEDirectForcingKinematics"),
        ib_method_ops,
        /*part*/ 0,
        /*register_for_restart*/ true);
    ib_method_ops->registerDirectForcingKinematics(dfk, /*part*/ 0);

    // 刚体自由度：先自由平移+绕 z 轴转动（与原 eel2d 对齐）
    IBTK::FreeRigidDOFVector solve_dofs;
    solve_dofs << 1, 1, 0, 0, 0, 1;
    dfk->setSolveRigidBodyVelocity(solve_dofs);

    // 4) 注册 eel 运动学函数（见 Eel2dKinematics.*）
    auto* eel_ctx = new Eel2dKinematicsCtx(/*input_db 子库*/);
    dfk->registerKinematicsFunction(&eel_com_kinematics_fcn, eel_ctx);

    // 5) 初始化 FE 方程系统
    ib_method_ops->initializeFEEquationSystems();

    // 6) 初始化 hierarchy，主时间推进循环
    //    在每一步里调用 eel_ctx->updateShapeAndMeshVelocity(t, ...)
    //    并写出可视化数据。
}
```

### 说明

- `registerDirectForcingKinematics()` 和 `setSolveRigidBodyVelocity()` 负责处理刚体自由度与约束求解。
- 变形部分在“最小版本”里先通过 `Eel2dKinematics` 维护一个目标形状场（见下节）。

---

## 2.2 `Eel2dKinematics.h`：上下文与接口

```cpp
#pragma once
#include <Eigen/Core>
#include <vector>

struct Eel2dKinematicsCtx
{
    // 参数：从 input2d 读取（可复用 ConstraintIB/eel2d 的表达式）
    double init_angle = 0.0;
    double wave_omega = 0.785 / 0.125;

    // 可选：muParser 对象，保存 body_shape_equation / deformation_velocity_function
    // mu::Parser body_shape_parser;
    // mu::Parser vel_parser_x, vel_parser_y;

    // 每个节点的参考弧长参数 s_i、目标位置 X_target、目标速度 U_target
    std::vector<double> s;
    std::vector<Eigen::Vector2d> X_target;
    std::vector<Eigen::Vector2d> U_target;

    void initializeFromInput(/*db*/);
    void initializeLagrangianParam(/*mesh / node map*/);
    void updateShapeAndMeshVelocity(double time,
                                    const Eigen::Vector3d& Xcom,
                                    const Eigen::Vector3d& Wcom);
};

// 只负责 COM 刚体速度（给 IBFEDirectForcingKinematics）
void eel_com_kinematics_fcn(double data_time,
                            Eigen::Vector3d& U_com,
                            Eigen::Vector3d& W_com,
                            void* ctx);
```

---

## 2.3 `Eel2dKinematics.cpp`：关键函数骨架

```cpp
#include "Eel2dKinematics.h"
#include <cmath>

void Eel2dKinematicsCtx::updateShapeAndMeshVelocity(double t,
                                                    const Eigen::Vector3d& Xcom,
                                                    const Eigen::Vector3d& Wcom)
{
    // 1) 计算 body frame 下每个节点目标位形（复用 eel2d 原公式思想）
    // y_base(s,t) = A(s) * sin(2*pi*s - omega*t)
    // U_def(s,t)  = d/dt[y_base] * n_hat

    for (size_t i = 0; i < s.size(); ++i)
    {
        const double si = s[i];
        const double A  = 0.125 * ((si + 0.03125) / 1.03125);
        const double phase = 2.0 * M_PI * si - (0.785 / 0.125) * t;

        const double y = A * std::sin(phase);
        const double dydt = -0.785 * ((si + 0.03125) / 1.03125) * std::cos(phase);

        // body frame 法向（最小版可先近似取 y 方向）
        Eigen::Vector2d n_hat(0.0, 1.0);

        Eigen::Vector2d Xb(si, y);
        Eigen::Vector2d Ub = dydt * n_hat;

        // 2) 叠加刚体旋转（由 Xcom, Wcom 给出）
        const double theta = /* 从姿态管理器读取 z 轴角，最小版可先 0 */ 0.0;
        Eigen::Matrix2d R;
        R << std::cos(theta), -std::sin(theta),
             std::sin(theta),  std::cos(theta);

        X_target[i] = R * Xb + Xcom.head<2>();

        // U = U_com + R*Ub + W x r
        Eigen::Vector2d r = R * Xb;
        Eigen::Vector2d Wxr(-Wcom[2] * r[1], Wcom[2] * r[0]);
        U_target[i] = U_com.head<2>() + R * Ub + Wxr;
    }

    // 3) 把 X_target/U_target 写入 FE 对应向量/系统（由你在 example.cpp 中接线）
}

void eel_com_kinematics_fcn(double /*data_time*/,
                            Eigen::Vector3d& U_com,
                            Eigen::Vector3d& W_com,
                            void* /*ctx*/)
{
    // 最小版：COM 速度由流体耦合自然求解 => 给 0，交给 DFK 的刚体求解器
    U_com.setZero();
    W_com.setZero();
}
```

---

## 3. `input2d` 最小参数块（新增/对齐）

```text
IBFEMethod {
  max_levels = MAX_LEVELS
  use_consistent_mass_matrix = TRUE
}

EelIBFEDirectForcingKinematics {
  rho_solid = RHO
  # 可按需打开输出
  print_output = TRUE
}

EelKinematics {
  initial_angle_body_axis_0       = 0.0
  body_shape_equation             = "0.125*((X_0+0.03125)/1.03125)*sin(2*PI*X_0-(0.785/0.125)*T)"
  deformation_velocity_function_0 = "(-0.785*((X_0+0.03125)/1.03125)*cos(2*PI*X_0-(0.785/0.125)*T))*N_0"
  deformation_velocity_function_1 = "(-0.785*((X_0+0.03125)/1.03125)*cos(2*PI*X_0-(0.785/0.125)*T))*N_1"
}
```

> 这样可以最大程度重用原 `ConstraintIB/eel2d` 的参数含义，便于对比。

---

## 4. 最小验证流程（务实版）

1. **几何验证**：先在无流体或弱耦合下输出 Lagrangian 点，看波形是否和原例子一致。
2. **动力学验证**：打开完整耦合，比较 3 个量：
   - COM 前进速度随时间；
   - 尾迹主涡尺度；
   - 一周期位移（stride length）。
3. **参数回归**：仅调 `dt`、网格、`rho_solid`，避免一次改太多。

---

## 5. 两个常见坑（提前规避）

- **坑 1：把 DFK 当成“自动给形变”**
  - `IBFEDirectForcingKinematics` 主要是刚体与约束接口；你仍需在外部维护鱼体目标形变场并写入 FE 向量（最小版可用惩罚/目标点策略）。

- **坑 2：上来就做主动应力**
  - 主动应力（PK1）是更物理，但最小复现不需要第一步就上。
  - 建议先完成“运动学复现”，再升级到“机理复现”。

---

## 6. 升级路径（从最小版到科研版）

- V1（本清单）：DFK + 目标位形/速度回放。
- V2：加入 FE 弹性本构，形变由“目标驱动 + 弹性”共同决定。
- V3：把目标驱动替换成时空主动应力（`registerPK1StressFunction`），形成自主推进机理模型。

