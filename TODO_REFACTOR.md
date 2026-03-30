# Lux Pure Rust Refactor TODO

本文档整理当前对 `luxpy` 的探索结论、纯 Rust 重构优先级、阶段目标与验收要求，作为后续执行参考。

## 1. 当前仓库状态

- Rust crate 当前已完成一版 `P0` 基础数值内核，已具备：
  - 波长网格生成 `getwlr`
  - 波长间距计算 `getwld`
  - 单谱数据模型 `Spectrum`
  - 批量谱数据模型 `SpectralMatrix`
  - 最小版线性 `cie_interp`
    - 线性插值
    - 线性外推
    - 可选负值裁剪
  - 最小版 `spd_normalize`
    - `max`
    - `area`
    - `lambda`
    - `ru`
    - `pu`
    - `qu`
  - 公开观察者 / CMF API
    - `xyzbar`
    - `vlbar`
  - 基础积分链路
    - `spd_to_power`
    - `spd_to_xyz`
    - `spd_to_ler`
  - 上述积分链路均支持单谱，`spd_to_xyz` / `spd_to_ler` 已支持批量谱
  - 嵌入式标准观察者：
    - `1931_2`
    - `1964_10`
  - CIE 191:2010 mesopic 支持：
    - `get_cie_mesopic_adaptation`
    - `vlbar_cie_mesopic`
  - 参考光源：
    - `blackbody`
    - `daylightlocus`
    - `daylightphase`
    - `cri_ref`（默认 CIE Ra 路径）
    - `standard_illuminant(name, wl_grid)` 首批 registry
  - CCT / Duv：
    - `xyz_to_cct`
    - `cct_to_xyz`
  - 常用颜色空间变换：
    - `XYZ <-> Yxy`
    - `XYZ <-> Yuv`
    - `XYZ <-> Lab`
    - `XYZ <-> Luv`
    - `XYZ <-> LMS`
    - `XYZ <-> sRGB`
  - 色适应 / 色差：
    - `cat_apply`
    - `cat_apply_mode`
    - `cat_apply_with_conditions`
    - `cat_apply_context`
    - `CatViewingConditions`
    - `CatContext`
    - `CatAdapter`
    - `deltaE`
- 当前相关文件：
  - `src/spectrum.rs`
  - `src/photometry.rs`
  - `src/color.rs`
  - `src/illuminants.rs`
  - `src/error.rs`
  - `data/cmfs/`
  - `data/spds/`
- 当前测试状态：
  - Rust 单测：130
  - Python parity 集成测试：1
  - `cargo test` 全通过

## 2. luxpy 顶层探索结论

`luxpy` 顶层入口在 `luxpy/luxpy/__init__.py`，实际暴露的是一个较完整的照明与色科学计算平台，不只是几个基础函数。

### 2.1 顶层核心计算域

- 光谱基础
  - `getwlr`
  - `getwld`
  - `spd_normalize`
  - `cie_interp`
  - `spd`
  - `SPD`
- 标准观察者 / CMF
  - `_CMF`
  - `xyzbar`
  - `vlbar`
  - `vlbar_cie_mesopic`
- SPD 积分链路
  - `spd_to_xyz`
  - `spd_to_ler`
  - `spd_to_power`
- 参考光源
  - `blackbody`
  - `daylightlocus`
  - `daylightphase`
  - `cri_ref`
- CCT / Duv
  - `xyz_to_cct`
  - `cct_to_xyz`
- 常用颜色空间变换
  - `XYZ <-> Yxy`
  - `XYZ <-> Yuv`
  - `XYZ <-> Lab`
  - `XYZ <-> Luv`
  - `XYZ <-> LMS`
  - `XYZ <-> sRGB`
- 色差
  - `deltaE`

### 2.2 色貌模型

`luxpy.color.cam` 是一个完整子系统，包含：

- CIECAM02
- CAM16
- CAM02-UCS / CAM16-UCS
- CAM15u
- CAM18sl
- ZCAM / Jabz
- SWW16
- 大量 `xyz_to_jab*` / `*_to_xyz` wrapper

结论：色貌模型必须晚于基础 XYZ、CAT、白点、观察条件建模。

### 2.3 色适应

`luxpy.color.cat` 支持多种 CAT：

- HPE
- CAT02
- CAT16
- Bradford
- Sharp
- CMCCAT2000
- Kries / Judd / Bianco 等

核心能力是对应色计算与适应度计算。

### 2.4 标准观察者与 CMF

`luxpy._CMF['types']` 当前包含：

- `1931_2`
- `1964_10`
- `2006_2`
- `2006_10`
- `2015_2`
- `2015_10`
- `1931_2_judd1951`
- `1931_2_juddvos1978`
- `1951_20_scotopic`
- `cie_std_dev_obs_f1`

每套不仅有 `bar`，还有：

- `K`
- `M`
- `N`

结论：Rust 端后续不能把观察者仅当作三列 CMF 数据，还要考虑附属常数与矩阵。

### 2.5 工具箱层能力

顶层还会挂载多个 toolbox 模块：

- `photbiochem`
- `indvcmf`
- `spdbuild`
- `hypspcim`
- `iolidfiles`
- `spectro`
- `rgb2spec`
- `dispcal`
- `sherbrooke_spectral_indices`
- `spectral_mismatch_and_uncertainty`

说明：

- `hypspcim` 对应高光谱图像模拟
- `photbiochem` 对应 alpha-opic、EDI、DER、ELR、蓝光危害、昼夜节律
- `indvcmf` 对应个体观察者模型
- `dispcal` / `spectro` / `iolidfiles` 更像外围工程工具链

## 3. 已识别的兼容性风险

### 3.1 `getwld()` 返回类型差异

`luxpy.getwld()`：

- 等间距波长时返回标量
- 非等间距时返回数组

当前 Rust `getwld()`：

- 一律返回 `Vec<f64>`

结论：需要尽早确定 Rust API 策略：

- 方案 A：内部统一表示，外部 API 不完全仿 LuxPy
- 方案 B：外部接口尽量兼容 LuxPy 语义

该决策会影响后续：

- 积分
- 插值
- 归一化
- 批量谱表示

### 3.2 单谱 / 多谱表示差异

`luxpy` 大量函数接受 `(N+1, M)` 形式矩阵：

- 第 0 行是波长
- 后续每一行是一条光谱

当前 Rust 偏向单条 `Spectrum`。

结论：P0 阶段必须设计统一的单谱与批量谱抽象。

### 3.3 插值语义复杂

`cie_interp()` 不只是简单线性插值，还包括：

- 按数据类型选择默认插值方法
  - `spd`
  - `cmf`
  - `rfl`
  - `none`
- 外推策略
- 是否允许负值
- 可选 Sprague / cubic / quadratic / linear
- 对数域插值 / 外推选项

结论：该函数应视为独立基础设施，不应作为零散 helper 实现。

## 4. 重构优先级矩阵

| 功能域 | 代表功能 | 依赖 | Python 对拍难度 | Rust 实现难度 | 优先级 | 说明 |
|---|---|---|---|---|---|---|
| 光谱底座 | `Spectrum/SPD`、`getwlr`、`getwld`、批量谱数据结构 | 无 | 低 | 中 | P0 | 所有后续能力的承载层 |
| 光谱插值与归一化 | `cie_interp`、`spd_normalize` | 光谱底座 | 中 | 高 | P0 | `luxpy` 核心语义之一 |
| 标准观察者与 CMF 数据 | `_CMF`、`xyzbar`、`vlbar`、`vlbar_cie_mesopic` | 光谱底座、插值 | 低 | 中 | P0 | 后续 CCT/CAT/CAM/CRI 都依赖 |
| 基础积分链路 | `spd_to_power`、`spd_to_ler`、`spd_to_xyz` | 上述全部 | 低到中 | 中 | P0 | 核心数值计算心脏 |
| 参考光源 | `blackbody`、`daylightlocus`、`daylightphase`、`cri_ref` | `spd_to_xyz`、CMF、插值 | 中 | 中到高 | P1 | 功能价值高，也是 CRI 基础 |
| 常用颜色空间变换 | `XYZ<->Yxy/Yuv/Lab/Luv/LMS/sRGB` | `spd_to_xyz`、白点/CMF | 低 | 中 | P1 | 使用频率高 |
| CCT / Duv | `xyz_to_cct`、`cct_to_xyz` | `spd_to_xyz`、颜色空间、参考光源 | 中到高 | 高 | P1 | 建议先做一条主算法 |
| 固定标准光源数据集 | CIE `A`、`D50/D55/D65/D75`、`F1..F12`、CIE LED 系列 | 光谱底座、插值、统一 illuminant API | 低到中 | 中 | P1.5 | 主要是数据资产、命名检索与重采样接口 |
| 色适应 | `cat.apply()`、适应度函数 | XYZ 变换、观察者矩阵 | 中 | 中到高 | P2 | 是 CAM 前置 |
| 色差 | `deltaE` | Lab/UCS 等颜色空间 | 低 | 中 | P2 | 依赖清晰，可中期交付 |
| 色貌模型 | CIECAM02、CAM16、CAM-UCS、ZCAM、CAM15u、CAM18sl | CAT、XYZ、观察条件 | 高 | 很高 | P3 | 体系大，不宜过早展开 |
| 显色评价 | `cri`、TM-30、`Rf/Rg`、CIE Ra | 参考光源、XYZ、色空间、CAT、数据库 | 高 | 很高 | P3 | 强依赖全链路稳定 |
| 光生物与节律 | `photbiochem`、alpha-opic、EDI、DER、ELR、BLH | SPD 积分、作用谱数据库 | 中 | 中 | P3 | 独立性较好，可后置并行 |
| 个体观察者模型 | `indvcmf` | CMF、矩阵、模型参数 | 高 | 高 | P4 | 研究型扩展，后置 |
| 高光谱图像 | `hypspcim` | 光谱底座、反射率库、颜色变换 | 高 | 很高 | P4 | 工程量大，后置 |
| 外围工具链 | `dispcal`、`rgb2spec`、`spectro`、`iolidfiles` | 各自独立 | 高 | 很高 | P4 | 更适合后续独立 crate / 工具层 |

## 5. 分阶段执行计划

### Phase P0: 基础数值内核

目标：形成可承载 `luxpy` 核心光谱数值语义的 Rust 底座。

计划：

- [x] 设计统一的单谱 / 多谱数据模型
- [x] 明确 `getwld` 的 Rust 语义
- [x] 重构波长网格与 spacing 表示
- [x] 实现最小版 `cie_interp`
- [x] 实现最小版 `spd_normalize`
- [x] 扩展基础观察者表示
- [x] 实现 `xyzbar`
- [x] 实现 `vlbar`
- [x] 实现 `vlbar_cie_mesopic`
- [x] 将 `spd_to_power` 重构到统一积分框架
- [x] 实现 `spd_to_ler`
- [x] 实现 `spd_to_xyz`

当前说明：

- 这里的 `cie_interp` 仍是“最小公共版”，尚未覆盖 LuxPy 的完整数据类型策略、Sprague、二次 / 三次外推、log 插值等全部语义。
- `spd_normalize` 当前也先覆盖最常用单谱模式，已经能满足基础内核验收，但尚未完全复刻 LuxPy 的所有列表参数和外围行为。
- 更全量 `_CMF` 集合以及更完整插值策略仍待补齐。
- 常用颜色空间当前已实现 `Yxy / Yuv / Lab / Luv / LMS / sRGB`。
- `Lab / Luv` 当前先采用“显式白点输入”的 Rust API，后续再决定是否补充默认参考白点与 illuminant 封装层。
- `LMS` 当前同时支持“显式 3x3 矩阵”与“基于 `Observer` 的默认矩阵”两层 API，后续新观察者接入时应沿用同一模式。
- `sRGB` 当前按 LuxPy 的 IEC:61966 语义实现，保持 `XYZ -> sRGB` 输出 `0..255` 浮点 RGB，后续如需 UI / image pipeline 可再补 `u8` 包装层。

P0 验收完成标志：

- Rust 端可稳定对拍：
  - [x] `getwlr`
  - [x] `getwld`
  - [x] 最小版 `cie_interp`
  - [x] 最小版 `spd_normalize`
  - [x] `xyzbar`
  - [x] `vlbar`
  - [x] `spd_to_power`
  - [x] `spd_to_ler`
  - [x] `spd_to_xyz`

### Phase P1: 参考光源 + CCT

计划：

- [x] 实现 `blackbody`
- [x] 实现 `daylightlocus`
- [x] 实现 `daylightphase`
- [x] 实现 `cri_ref`
- [x] 实现 `XYZ <-> Yxy`
- [x] 实现 `XYZ <-> Yuv`
- [x] 实现 `XYZ <-> Lab`
- [x] 实现 `XYZ <-> Luv`
- [x] 实现 `XYZ <-> LMS`
- [x] 实现 `XYZ <-> sRGB`
- [x] 选择一条主路径实现 `xyz_to_cct`
- [x] 实现 `cct_to_xyz`

说明：

- `xyz_to_cct` 不要一开始追平全部 Robertson / Ohno / Li / Zhang 分支
- 先做一条稳定、可验收、可维护的主算法
- 颜色空间转换层已开始收拢为统一 API，后续新增 `LMS / sRGB` 时应复用相同的显式参数风格或一次性升级默认白点策略，避免接口分裂。
- P1 当前状态可以视为“颜色空间子集、默认参考光源链路与 CCT 主链路已完成”。

### Phase P1.5: 固定标准光源数据集

目标：补齐依赖固定 SPD 表的标准光源数据，并统一到可重采样的 illuminant API。

计划：

- [x] 设计标准光源命名与查找接口
  - [x] 统一 `standard_illuminant(name, wl_grid)` 风格入口
  - [ ] 明确命名规范与别名收口
  - [x] 提供基础错误处理
- [x] 接入固定标准光源数据集首批集合
  - [x] CIE `A`
  - [x] CIE `D50`
  - [x] CIE `D55`
  - [x] CIE `D65`
  - [x] CIE `D75`
  - [x] CIE `F1..F12`
  - [x] CIE LED 系列
- [x] 统一重采样语义
- [x] 为首批固定标准光源建立 Python 对拍基线

说明：

- 该阶段本质上是“数据资产 + API 收口”，不是新的核心数值算法。
- 实现时不要继续堆散装函数，优先做统一 registry / lookup 层。
- `blackbody/daylightphase/cri_ref` 属于“生成型参考光源”，本阶段属于“表驱动固定光源”，两类入口要并存但 API 风格应一致。
- 当前状态：首批 registry 已完成，相关 CSV 已迁入仓库自有 `data/spds/`，不再依赖 `luxpy/` 目录下的数据文件。
- 当前未完成项主要是别名体系、更多标准光源数据集，以及后续非默认 CRI / TM-30 用到的混合参考光源策略。

### Phase P2: 色适应 + 色差

计划：

- [x] 实现 `cat.apply()` 主路径
- [x] 实现适应度计算
- [x] 实现基础 CAT 变换族
  - [x] Bradford
  - [x] CAT02
  - [x] CAT16
  - [x] Sharp
  - [x] Bianco
  - [x] CMC
  - [x] Kries
  - [x] Judd1945
  - [x] Judd1945Cie016
  - [x] Judd1935
- [x] 实现 `deltaE`

说明：

- 当前 `CAT` 已提供一步 von Kries 主路径，公开 API 以 `XYZ + source white + target white + transform + D` 形式暴露。
- 当前已实现的变换矩阵包括 `Bradford`、`CAT02`、`CAT16`、`Sharp`、`Bianco`、`CMC`、`Kries`、`Judd1945`、`Judd1945Cie016`、`Judd1935`。
- 当前已补上基于环境参数的适应度 `D` 计算，以及基于观察条件的上层适配入口。
- 当前已补上 `CatMode` 策略层，覆盖 `1>2`、`1>0`、`0>2`、`1>0>2`。
- 当前已补上更高层观察条件工具，包括 `CatViewingConditions`、`CatContext`，以及 `cat_apply_context`。
- 当前已补上可复用的批量/矩阵风格 CAT 封装，包括预编译 `CatAdapter`、`cat_compile*`、`Tristimulus::cat_apply_adapter()`、`TristimulusSet::cat_apply_adapter()`。
- 当前尚未实现的是更多文献模型，以及更完整的高层 CAT utility 生态（例如更系统的预计算/缓存策略与进一步的批量工作流封装）。
- 当前已落地 `deltaE` 首批主路径：
  - [x] `CIE76`
  - [x] `CIEDE2000`
- 当前 Rust API 以 `XYZ + 白点` 作为公开入口，内部再转换到 `Lab` 执行公式计算；后续如需贴近 `luxpy` 的更泛化 `DE_cspace` 族接口，再在此基础上扩展。

### Phase P3: CAM + CRI + photobiochem

计划：

- [x] 实现 CIECAM02 首批前向/反向主路径
- [x] 实现 CAM16 首批前向/反向主路径
- [x] 实现 CAM02-UCS / CAM16-UCS 首批前向/反向主路径
- [x] 实现部分 `xyz_to_jab*` wrapper
- [x] 实现 CIE Ra
- [x] 实现 CIE Rf / Rg
- [x] 实现 IES TM-30 风格别名入口
- [x] 实现首批 TM-30 结果对象
- [ ] 实现 `photbiochem` 基础能力
  - [ ] alpha-opic irradiance
  - [ ] EDI
  - [ ] DER
  - [ ] ELR
  - [ ] blue light hazard

说明：

- 当前已先落地 CAM 前置能力，包括 `cam_naka_rushton`、`CamModel`、`CamSurround`、`CamViewingConditions`，以及 `cam16_viewing_conditions()` / `ciecam02_viewing_conditions()`。
- 当前这些能力主要负责观察条件、白点适应与响应压缩等底层公共计算，为后续 `CIECAM02` / `CAM16` 正式前向模型提供共享基础。
- 当前已开始接入首批正式前向模型，提供 `cam_forward()`、`cam16_forward()`、`ciecam02_forward()` 以及 `CamAppearance`，可输出 `J/Q/C/M/s/h` 与 `aM/bM`、`aC/bC`。
- 当前已补上 `J+aM+bM -> XYZ` 反向路径，以及 `CAM02-UCS / CAM16-UCS` 的 `J'a'b' <-> XYZ` 主路径。
- 当前已补上 `Tristimulus/TristimulusSet` 风格 CAM / CAM-UCS wrapper。
- 当前已补上更高层 `xyz_to_jab*` / `jab*_to_xyz` 便捷 API，以及统一 `CamCoordinates / CamSpace` 坐标入口。
- 当前已补上 `CIE Ra`、`CIE Rf / Rg` 主路径，样品数据已迁入仓库自有 `data/rfls/`，不再依赖本地 `luxpy` 安装路径。
- 当前已补上 `IES TM-30` 风格别名入口，直接复用 `CIE 224` 的 `Rf / Rg` 数值主路径。
- 当前已补上首批 `TM-30` 结果对象，包括 `Rf/Rg/Rfi`、`DEi`、8 个 hue bin 的平均 `Jab` 点，以及 bin-level 的 local fidelity、chroma shift、hue shift。
- 当前尚未实现的是更完整的属性输入组合反解，以及 `photbiochem` 基础能力。

### Phase P4: 高级扩展与工具箱

计划：

- [ ] `indvcmf`
- [ ] `hypspcim`
- [ ] `dispcal`
- [ ] `rgb2spec`
- [ ] `spectro`
- [ ] `iolidfiles`
- [ ] `spectral_mismatch_and_uncertainty`

## 6. Python 对拍验收要求

以后每个功能都必须按“双轨验收”执行：

- Rust 单元测试 / 集成测试
- Python `luxpy` 对拍

### 6.1 对拍原则

- 每个功能至少覆盖：
  - 固定教科书样例
  - 边界样例
  - 批量 / 矩阵样例
- 不只比单个输出值，要比一组固定 case
- 对拍脚本统一使用：
  - `luxpy/.venv/bin/python`
- 对拍环境统一设置：
  - `MPLCONFIGDIR=/tmp/mpl`

### 6.2 建议目录

- `tests/python_ref/`
  - Python 参考计算脚本
- `tests/parity/`
  - Rust 对拍测试

### 6.3 误差建议

- 基础积分 / 颜色空间：
  - 目标误差优先控制在 `1e-9` 到 `1e-6`
- CCT / CAM / CRI：
  - 根据算法稳定性单独设阈值

### 6.4 首批必须建立的 Python 参考基线

- [x] `getwlr`
- [x] `getwld`
- [x] `spd_to_power(ru/pu/qu)`
- [x] `spd_to_xyz`
- [x] `spd_to_ler`
- [x] `xyzbar`
- [x] `vlbar`
- [x] 最小版 `cie_interp`
- [x] 最小版 `spd_normalize`
- [x] `blackbody`
- [x] `daylightphase`
- [x] `cri_ref`
- [x] `xyz_to_cct`
- [x] `cct_to_xyz`
- [x] 固定标准光源数据集首批 registry
  - [x] CIE `A`
  - [x] CIE `D50/D55/D65/D75`
  - [x] CIE `F1..F12`
  - [x] CIE LED 系列

## 7. 已确认的 Python 参考值

以下结果已通过本地 `luxpy` 虚拟环境确认，可作为初始 sanity check：

- `getwld([400,410,420]) = 10.0`
- `getwld([400,410,430]) = [10.0, 15.0, 20.0]`
- `spd_to_power([[400,410,420],[1,2,3]], 'ru') = 60.0`
- `spd_to_power([[555,556],[1,1]], 'pu', '1931_2') = 1365.9061258134`
- `spd_to_power([[500,510],[1,1]], 'qu') = 5.084457733218137e19`
- `spd_to_ler([[555,556],[1,1]], '1931_2') = 682.9530629067`
- `spd_to_ler([[555,556],[1,1]], '1964_10') = 683.2161845600001`
- `spd_to_xyz([[555,556],[1,1]], '1931_2') = [52.02102730660652, 100.0, 0.5527195523559263]`
- `spd_to_xyz([[555,556],[1,1]], '1931_2', relative=False) = [710.558398692, 1365.9061258134, 7.549630224198]`
- `spd_to_xyz([[555,556],[1,1]], '1964_10') = [62.535069638997825, 100.0, 0.09015048427119186]`
- `xyzbar('1931_2')` 与 `vlbar('1931_2')` 已建立 555 nm 和线性重采样对拍基线
- `xyzbar('1964_10')` 与 `vlbar('1964_10')` 已建立 555 nm 对拍基线
- `xyz_to_cct([100,100,100]) ≈ (5455.5 K, -0.0044233)`
- `blackbody(6500, wl=[360,365,1])` 已可生成参考谱
- `daylightphase(6500, wl=[360,365,1])` 已可生成参考谱
- `cri_ref([3000,6500], wl=[360,365,1])` 已可生成参考谱
- `cct_to_xyz(6500)` 已建立对拍基线
- `standard_illuminant('A'/'D50'/'D65'/'F4'/'LED_B1', ...)` 已建立对拍基线
- 标准光源 CSV 当前位于仓库 `data/spds/`

## 8. 建议的近期执行顺序

近期不应同时铺开太多方向，建议严格按下面顺序推进：

1. [x] 重构 `Spectrum/SPD` 与批量谱抽象
2. [x] 实现最小版 `cie_interp`
3. [x] 实现最小版 `spd_normalize`
4. [x] 扩展 `_CMF` / `xyzbar` / `vlbar`
5. [x] 打通 `spd_to_xyz`
6. [x] 建立当前阶段 Python 对拍基线
7. [x] 进入 `blackbody/daylightphase/cri_ref`
8. [x] 再进入 `xyz_to_cct/cct_to_xyz`
9. [x] 然后进入固定标准光源数据集层首批 registry
10. [ ] 再进入更高层的 `CAT` 工具层

当前建议：

- `blackbody / daylightphase / daylightlocus / cri_ref` 默认链路已完成
- `xyz_to_cct / cct_to_xyz` 主路径已完成
- 固定标准光源首批 registry 已完成
- 光谱静态数据已经开始沉淀到仓库自有 `data/` 目录，后续新增标准光源应沿用同一组织方式
- `deltaE` 首批主路径已完成
- `CAT` 的一步主路径、`Bradford/CAT02/CAT16/Sharp/Bianco/CMC/Kries/Judd` 族、`D` / 观察条件入口、模式层、`CatContext`、`CatAdapter` 与 `cat_compile*` 已完成
- 下一阶段优先进入 `photbiochem` 基础能力
- 固定标准光源后续扩展仍应沿用统一 illuminant registry，而不是分散新增 `f1()`、`led_b1()` 一类 API

## 9. 暂不做的事项

在基础内核稳定前，暂不优先处理：

- 全量 CAM 族
- 全量 CRI / TM-30 图形输出
- 高光谱图像模拟
- 仪器接口
- 显示器校准
- 立体显示 / VR viewer

这些都依赖更底层的正确性，过早进入会明显放大验证成本。
