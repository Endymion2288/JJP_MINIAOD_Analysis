# JJP SPS MiniAOD Gen-Level Analysis

分析 MiniAOD 中 Jpsi1, Jpsi2, phi 粒子之间的 delta eta 和 delta phi 关联（SPS 选择标准）。

## 分析目标
- 从 MiniAOD 的 `prunedGenParticles` 集合中提取 J/psi (pdgId=443) 和 phi (pdgId=333) 粒子
- 根据 SPS（Single Parton Scattering）物理要求选取粒子：
  - **Jpsi1, Jpsi2, phi** 三者必须来源于**同一个部分子对撞过程**（共享共同祖先）
  - **Jpsi1 和 Jpsi2** 按 pT 排序：Jpsi1 pT > Jpsi2 pT
  - **phi**：选取来源于 gluon 衰变的 phi（gluon 来自部分子对撞产生）
  - 若一个 event 中存在多个候选组合，选择**总 pT（Jpsi1 + Jpsi2 + phi）最大**的组合
- 计算并绘制：
  1. Delta rapidity 和 delta phi 的 1D 直方图
  2. Delta rapidity vs delta phi 的 2D 关系图
  3. 运动学分布（pT, eta）

## 与 DPS_2 的区别
| 选择标准 | DPS_2 | SPS |
|---------|-------|-----|
| Jpsi1-Jpsi2 关系 | 共享共同祖先 | 三者共享共同祖先 |
| phi 选择 | 最高 pT 的来自 gluon 的 phi | 最高 pT 的来自 gluon 的 phi，且与两个 Jpsi 共享共同祖先 |
| 候选组合选择 | Jpsi pair 总 pT 最大 | 三者总 pT 最大 |

## 文件说明
- `analyze_gen_correlations.py`: 主分析脚本（FWLite-based Python 分析）
- `plot_results.py`: ROOT 输出文件的绘图脚本
- `run_analysis.sh`: 一键运行分析和绘图的脚本

## 使用方法

```bash
cd /eos/home-x/xcheng/learn_MC/loopmix_pythia/CMSSW_14_0_18/src
eval `scramv1 runtime -sh`  # 或 cmsenv
cd JJP_SPS_MINIAOD_Analysis

# 运行分析（可调整参数）
python3 analyze_gen_correlations.py \
    -i /eos/user/x/xcheng/learn_MC/JJP_SPS_MC_output/MINIAOD/ \
    -o gen_correlation_histograms_sps.root \
    -n -1 \
    --max-files -1 \
    -j 4  # 并行 worker 数

# 生成绘图
python3 plot_results.py \
    -i gen_correlation_histograms_sps.root \
    -o plots
```

### 命令行参数

**analyze_gen_correlations.py**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i, --input-dir` | MiniAOD 输入目录 | `/eos/user/x/xcheng/learn_MC/JJP_SPS_MC_output/MINIAOD/` |
| `-o, --output` | 输出 ROOT 文件名 | `gen_correlation_histograms_sps.root` |
| `-n, --max-events` | 最大处理事件数（-1 表示全部） | -1 |
| `--max-files` | 最大处理文件数（-1 表示全部） | -1 |
| `-j, --jobs` | 并行 worker 数 | 1 |

**plot_results.py**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i, --input` | 输入 ROOT 文件 | `gen_correlation_histograms_sps.root` |
| `-o, --output-dir` | 输出图像目录 | `plots` |

## 选择算法详解

### SPS 选择流程
1. 从 `prunedGenParticles` 提取所有 J/psi 和 phi 粒子（使用 `isLastCopy()` 或 status==2）
2. 按 pT 降序排序
3. 筛选来自 gluon 衰变的 phi 候选
4. 遍历所有 (Jpsi_i, Jpsi_j, phi_k) 组合：
   - 检查三者是否共享共同祖先（同一部分子对撞）
   - 若满足条件，计算总 pT = pT(Jpsi_i) + pT(Jpsi_j) + pT(phi_k)
5. 选择总 pT 最大的组合，其中 pT 较大的 Jpsi 为 Jpsi1

## 输出
- `gen_correlation_histograms_sps.root`: 包含所有直方图的 ROOT 文件
  - `h_dy_12`, `h_dy_1phi`, `h_dy_2phi`: Delta rapidity 1D 直方图
  - `h_dphi_12`, `h_dphi_1phi`, `h_dphi_2phi`: Delta phi 1D 直方图
  - `h2_dy_dphi_12`, `h2_dy_dphi_1phi`, `h2_dy_dphi_2phi`: 2D 关联图
  - `h_jpsi1_pt`, `h_jpsi2_pt`, `h_phi_pt`: pT 分布
  - `h_jpsi1_eta`, `h_jpsi2_eta`, `h_phi_eta`: eta 分布
- `plots/`: 输出的 PDF/PNG 图像文件

## MiniAOD 输入文件
默认位于: `/eos/user/x/xcheng/learn_MC/JJP_SPS_MC_output/MINIAOD/`

## 依赖
- CMSSW_14_0_18（或兼容版本）
- FWLite
- ROOT with PyROOT
