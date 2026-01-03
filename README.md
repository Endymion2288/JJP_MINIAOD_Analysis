# JJP TPS MiniAOD Gen-Level Analysis

分析 TPS (Triple Parton Scattering) 样本 MiniAOD 中 Jpsi1, Jpsi2, phi 粒子之间的 delta y 和 delta phi 关联。

## 分析目标

从 MiniAOD 的 `prunedGenParticles` 集合中提取 J/ψ (pdgId=443) 和 φ (pdgId=333) 粒子，研究它们之间的运动学关联。

## TPS 粒子选取标准

- **Jpsi1**: pT 最高的 J/ψ
- **Jpsi2**: pT 次高的 J/ψ  
- **phi**: pT 最高的且来源于 gluon 衰变的 φ 介子

## 输出内容

### 1D 直方图
- `h_dy_12`: Jpsi1-Jpsi2 之间的 |Δy|
- `h_dy_1phi`: Jpsi1-Phi 之间的 |Δy|
- `h_dy_2phi`: Jpsi2-Phi 之间的 |Δy|
- `h_dphi_12`: Jpsi1-Jpsi2 之间的 |Δφ|
- `h_dphi_1phi`: Jpsi1-Phi 之间的 |Δφ|
- `h_dphi_2phi`: Jpsi2-Phi 之间的 |Δφ|
- `h_jpsi1_pt`, `h_jpsi2_pt`, `h_phi_pt`: pT 分布
- `h_jpsi1_eta`, `h_jpsi2_eta`, `h_phi_eta`: η 分布

### 2D 直方图
- `h2_dy_dphi_12`: Jpsi1-Jpsi2 的 Δy vs Δφ
- `h2_dy_dphi_1phi`: Jpsi1-Phi 的 Δy vs Δφ
- `h2_dy_dphi_2phi`: Jpsi2-Phi 的 Δy vs Δφ

## 文件说明

| 文件 | 描述 |
|------|------|
| `analyze_gen_correlations.py` | 主分析脚本，基于 FWLite 读取 MiniAOD |
| `plot_results.py` | 绘图脚本，生成 PDF/PNG 图像 |
| `run_analysis.sh` | 一键运行脚本 |

## 使用方法

### 方法 1: 使用一键脚本

```bash
cd /afs/cern.ch/user/x/xcheng/cernbox/learn_MC/loopmix_pythia/CMSSW_14_0_18/src/JJP_TPS_MINIAOD_Analysis
source run_analysis.sh
```

### 方法 2: 手动运行

```bash
# 设置 CMSSW 环境
cd /afs/cern.ch/user/x/xcheng/cernbox/learn_MC/loopmix_pythia/CMSSW_14_0_18/src
cmsenv

# 进入分析目录
cd JJP_TPS_MINIAOD_Analysis

# 运行分析 (可选参数)
python3 analyze_gen_correlations.py \
    -i /eos/user/x/xcheng/learn_MC/JJP_TPS_MC_output/MINIAOD/ \
    -o gen_correlation_histograms.root \
    -n -1 \
    --max-files -1 \
    -j 4  # 并行处理

# 绘图
python3 plot_results.py -i gen_correlation_histograms.root -o plots
```

## 命令行参数

### analyze_gen_correlations.py

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input-dir` | `/eos/user/x/xcheng/learn_MC/JJP_TPS_MC_output/MINIAOD/` | 输入 MiniAOD 文件目录 |
| `-o, --output` | `gen_correlation_histograms.root` | 输出 ROOT 文件 |
| `-n, --max-events` | -1 | 最大处理事件数 (-1 表示全部) |
| `--max-files` | -1 | 最大处理文件数 (-1 表示全部) |
| `-j, --jobs` | 1 | 并行处理的 worker 数量 |

### plot_results.py

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-i, --input` | `gen_correlation_histograms.root` | 输入直方图文件 |
| `-o, --output-dir` | `plots` | 输出图像目录 |

## MiniAOD 输入文件

位于: `/eos/user/x/xcheng/learn_MC/JJP_TPS_MC_output/MINIAOD/`

## 与 DPS 分析的区别

| 方面 | DPS | TPS |
|------|-----|-----|
| Jpsi1 选取 | 与 phi 的 gluon 祖先共享的 J/ψ | pT 最高的 J/ψ |
| Jpsi2 选取 | 剩余 pT 最高的 J/ψ | pT 次高的 J/ψ |
| phi 选取 | pT > 4 GeV 且来自 gluon 衰变，gluon 与某 J/ψ 共享祖先 | pT 最高的来自 gluon 衰变的 φ |
