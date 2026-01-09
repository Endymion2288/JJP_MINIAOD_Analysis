# JJP MiniAOD Gen-Level Analysis

分析 MiniAOD 中 J/ψ₁, J/ψ₂, φ 粒子之间的 delta eta 和 delta phi 关联。

## 分析目标

从 MiniAOD 的 `prunedGenParticles` 集合中提取以下粒子:
- **J/ψ₁** (pdgId = 443) - 按 pT 排序的第一个 J/ψ
- **J/ψ₂** (pdgId = 443) - 按 pT 排序的第二个 J/ψ
- **φ** (pdgId = 333)

根据不同的物理过程选择对应的粒子组合。

## 选择标准

| 模式 | 描述 |
|------|------|
| **SPS** | J/ψ₁, J/ψ₂, φ 来源于同一个部分子过程，φ 来自该过程产生的 gluon 衰变 |
| **DPS_1** | J/ψ₁ 和 φ 来源于同一个部分子过程（φ 来自 gluon），J/ψ₂ 来自另一个部分子过程 |
| **DPS_2** | J/ψ₁ 和 J/ψ₂ 来源于同一个部分子过程，φ 来自另一个部分子过程的 gluon 衰变 |
| **TPS** | 仅要求 φ 来自 gluon 衰变（不再要求三者来自不同部分子过程），候选按 pT 最高优先 |

## 输出内容

### 直方图
1. **1D 直方图**
   - Delta rapidity (|Δy|): J/ψ₁-J/ψ₂, J/ψ₁-φ, J/ψ₂-φ
   - Delta phi (|Δφ|): J/ψ₁-J/ψ₂, J/ψ₁-φ, J/ψ₂-φ
   - pT 分布: J/ψ₁, J/ψ₂, φ
   - η 分布: J/ψ₁, J/ψ₂, φ
   - rapidity (y) 分布: J/ψ₁, J/ψ₂, φ
   - 不变质量: M(J/ψ₁+J/ψ₂), M(J/ψ₁+φ), M(J/ψ₂+φ), M(J/ψ₁+J/ψ₂+φ)

2. **2D 直方图**
   - Δy vs Δφ: J/ψ₁-J/ψ₂, J/ψ₁-φ, J/ψ₂-φ

## 文件说明

| 文件 | 说明 |
|------|------|
| `analyze_gen_correlations.py` | 主分析脚本（FWLite-based Python 分析） |
| `plot_results.py` | ROOT 输出文件的绘图脚本 |
| `run_analysis.sh` | 单一模式分析运行脚本 |
| `run_all_modes.sh` | 运行所有模式的脚本 |

## 使用方法

### 环境设置

```bash
cd /eos/user/x/xcheng/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`  # 或 cmsenv
export X509_USER_PROXY=/afs/cern.ch/user/x/xcheng/x509up_u180107
cd JJPMCAnalyzer
```

### 运行单一模式分析

```bash
# 使用脚本运行（推荐，xrootd 输入）
./run_analysis.sh -m DPS_1 --max-files 10 -n 10000 -j 4

# 使用本地/EOS 输入目录（示例路径）
./run_analysis.sh -m DPS_1 -i /eos/user/x/xcheng/learn_MC/JJP_DPS_MC_output/MINIAOD/ --max-files 10 -n 10000 -j 4

# 或直接运行 Python 脚本
python3 analyze_gen_correlations.py \
    --xrootd-base /eos/ihep/cms/store/user/xcheng/MC_Production/output/JJP_DPS1 \
    -o gen_correlation_DPS_1.root \
    -m DPS_1 \
    -n -1 \
    --max-files -1 \
    -j 4
```

### 运行所有模式

```bash
./run_all_modes.sh
```

### 生成绘图

```bash
python3 plot_results.py \
    -i gen_correlation_DPS_1.root \
    -o plots_DPS_1 \
    -m DPS_1
```

## 命令行参数

### analyze_gen_correlations.py

| 参数 | 说明 |
|------|------|
| `-i, --input-dir` | 本地 MiniAOD 输入目录 |
| `--xrootd-base` | xrootd 上的输入目录基路径 |
| `--file-list` | 包含输入文件列表的文本文件 |
| `-o, --output` | 输出 ROOT 文件名 |
| `-m, --mode` | 选择模式：SPS, DPS_1, DPS_2, TPS |
| `-n, --max-events` | 最大处理事件数（-1 表示全部） |
| `--max-files` | 最大处理文件数（-1 表示全部） |
| `-j, --jobs` | 并行 worker 数（默认 1，顺序处理） |
| `--xrootd-server` | xrootd 服务器地址（默认 root://cceos.ihep.ac.cn） |

### run_analysis.sh

| 参数 | 说明 |
|------|------|
| `-m, --mode` | 选择模式（默认 DPS_1） |
| `--max-files` | 最大处理文件数 |
| `-n, --max-events` | 最大处理事件数 |
| `-j, --jobs` | 并行 worker 数 |
| `-o, --output-dir` | 输出目录（默认 output） |
| `--xrootd-base` | 外部指定 xrootd 输入路径 |
| `-i, --input-dir` | 本地或 EOS 输入目录（优先级高于 xrootd） |

## 输入文件位置

MiniAOD 文件位于 IHEP xrootd:
```
root://cceos.ihep.ac.cn//eos/ihep/cms/store/user/xcheng/MC_Production/output/
├── JJP_SPS/
├── JJP_DPS1/
│   ├── 0/output_MINIAOD.root
│   ├── 1/output_MINIAOD.root
│   └── ...
├── JJP_DPS2/
└── JJP_TPS/
```

## 输出示例

运行后会生成:
```
output/
├── gen_correlation_SPS.root
├── gen_correlation_DPS_1.root
├── gen_correlation_DPS_2.root
├── gen_correlation_TPS.root
├── plots_SPS/
│   ├── delta_y_comparison_SPS.pdf
│   ├── delta_phi_comparison_SPS.pdf
│   ├── correlation_2d_all_SPS.pdf
│   ├── pt_distributions_SPS.pdf
│   ├── eta_distributions_SPS.pdf
│   ├── rapidity_distributions_SPS.pdf
│   └── invariant_mass_SPS.pdf
├── plots_DPS_1/
├── plots_DPS_2/
└── plots_TPS/
```

## 粒子选择逻辑

### 共同祖先判定
通过追溯粒子的母粒子链，判断两个或多个粒子是否共享共同祖先，从而确定它们是否来自同一个部分子相互作用。

### φ 来自 gluon 衰变
遍历 φ 的母粒子链，检查是否存在 gluon (pdgId = 21) 祖先。

### J/ψ 选择
由于需要两个 J/ψ 粒子，程序要求事件中至少有 2 个 J/ψ。按 pT 从高到低排序后：
- J/ψ₁：pT 最高的 J/ψ
- J/ψ₂：pT 第二高的 J/ψ

### 选择优先级
对于每种模式，按 pT 从高到低排序候选粒子，选择满足条件的第一个组合。

## 注意事项

1. 需要有效的 X509 代理来访问 IHEP xrootd
2. 并行处理 (`-j > 1`) 使用 multiprocessing spawn 模式，避免 ROOT 在 fork 进程中的问题
3. 如果遇到文件读取错误，脚本会跳过该文件并继续处理
4. 与 JUP 分析不同，JJP 需要至少 2 个 J/ψ 粒子而非 1 个 J/ψ + 1 个 Υ
