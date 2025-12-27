# JJP DPS MiniAOD Gen-Level Analysis

分析 MiniAOD 中 Jpsi1, Jpsi2, phi 粒子之间的 delta eta 和 delta phi 关联。

## 分析目标
- 从 MiniAOD 的 `prunedGenParticles` 集合中提取 J/psi (pdgId=443) 和 phi (pdgId=333) 粒子
- 根据母粒子关系选取 Jpsi2 和 phi，要求它们来源于同一个部分子对撞过程
  - phi 来源于 gluon 衰变
  - 该 gluon 和 Jpsi2 在同一个部分子对撞中产生
- 计算并绘制：
  1. Delta eta 和 delta phi 的 1D 直方图
  2. Delta eta vs delta phi 的 2D 关系图

## 文件说明
- `analyze_gen_particles_cfg.py`: CMSSW EDAnalyzer 配置文件
- `JJP_GenAnalyzer.cc`: C++ EDAnalyzer 源代码
- `plot_results.py`: ROOT 输出文件的绘图脚本

## 使用方法

```bash
cd /eos/home-x/xcheng/learn_MC/loopmix_pythia/CMSSW_14_0_18/src
cmsenv
scram b -j 4
cmsRun JJP_DPS_MINIAOD_Analysis/analyze_gen_particles_cfg.py
python3 JJP_DPS_MINIAOD_Analysis/plot_results.py
```

## MiniAOD 输入文件
位于: `/eos/user/x/xcheng/learn_MC/JJP_DPS_MC_output/MINIAOD/`
