# JJP DPS MiniAOD Gen-Level Analysis (DPS_2)

分析 MiniAOD 中 Jpsi1, Jpsi2, phi 粒子之间的 delta eta 和 delta phi 关联。

## 分析目标
- 从 MiniAOD 的 `prunedGenParticles` 集合中提取 J/psi (pdgId=443) 和 phi (pdgId=333) 粒子
- 根据母粒子关系选取粒子：
  - **Jpsi1 和 Jpsi2**：来源于同一个部分子对撞过程（共享共同祖先）
    - 按 pT 排序：Jpsi1 pT > Jpsi2 pT
    - 若存在多个候选 Jpsi pair，选择总 pT 最大的组合
  - **phi**：选取 pT 最高的、且来源于 gluon 衰变的 phi
- 计算并绘制：
  1. Delta rapidity 和 delta phi 的 1D 直方图
  2. Delta rapidity vs delta phi 的 2D 关系图
  3. 运动学分布（pT, eta）

## 文件说明
- `analyze_gen_correlations.py`: 主分析脚本（FWLite-based Python 分析）
- `plot_results.py`: ROOT 输出文件的绘图脚本
- `run_analysis.sh`: 一键运行分析和绘图的脚本

## 使用方法

```bash
cd /eos/home-x/xcheng/learn_MC/loopmix_pythia/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`  # 或 cmsenv
cd JJP_DPS_MINIAOD_Analysis

# 运行分析（可调整参数）
python3 analyze_gen_correlations.py \
    -i /eos/user/x/xcheng/learn_MC/JJP_DPS_MC_output/MINIAOD/ \
    -o gen_correlation_histograms.root \
    -n -1 \
    --max-files -1 \
    -j 4  # 并行 worker 数

# 生成绘图
python3 plot_results.py \
    -i gen_correlation_histograms.root \
    -o plots
```

### 命令行参数
| 参数 | 说明 |
|------|------|
| `-i, --input-dir` | MiniAOD 输入目录 |
| `-o, --output` | 输出 ROOT 文件名 |
| `-n, --max-events` | 最大处理事件数（-1 表示全部） |
| `--max-files` | 最大处理文件数（-1 表示全部） |
| `-j, --jobs` | 并行 worker 数（默认 1，顺序处理） |

## MiniAOD 输入文件
位于: `/eos/user/x/xcheng/learn_MC/JJP_DPS_MC_output/MINIAOD/`

## 输出
- `gen_correlation_histograms.root`: 包含所有直方图的 ROOT 文件
- `plots/`: 输出的 PDF/PNG 图像文件
