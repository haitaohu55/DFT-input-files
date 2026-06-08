# V2O5 的 `g_cvq` 与 `A_ind` 计算教程

这份说明对应当前目录：

```text
/Users/huhaitao/Desktop/code/Material/19.V2O5
```

当前分工是：

```text
7_h5/      只放生成 h5 的 VASP 输入和输出
8_cvq/     只放从 vaspelph.h5 提取 g_cvq、计算 A_ind 的后处理脚本和结果
```

## 1. 目标量

我们要从 `vaspelph.h5` 中提取固定跃迁

```text
VBM: k_v = (0.4651162791, 0.5, 0.0)
CBM: k_c = (0.0, 0.0, 0.0)
q = k_c - k_v = (-0.4651162791, -0.5, 0.0)
  = (0.5348837209, 0.5, 0.0)  [fold 到 0-1 区间]
```

对应的电子-声子矩阵元为：

```text
g_cv,nu(q) = < psi_c,kc | dV_KS / dQ_nu,q | psi_v,kv >
```

其中 `nu = 1 ... 3N` 是固定 `q` 点上的声子支。当前 V2O5 原胞有 14 个原子，所以共有 42 支声子。

后续计算的系数是：

```text
A_ind = 5.3464012e-4 * (1 + 2 N_ph) * D_ph^2 / epsilon_ph * F0^(5/2)
```

多声子支时，正式总值按所有声子支求和：

```text
A_ind,total = 5.3464012e-4 * F0^(5/2)
              * sum_nu [(1 + 2 N_nu) * D_nu^2 / epsilon_nu]
```

不要只取最大 `|g|` 的那一支作为正式结果。可以把最大支当作近似或主导贡献分析，但总系数应全支求和。

## 2. 前置文件检查

先确认 `7_h5` 中已经有 VASP 生成的 h5 文件：

```bash
cd /Users/huhaitao/Desktop/code/Material/19.V2O5
ls -lh 7_h5/OUTPUT/vaspelph.h5 7_h5/OUTPUT/vaspout.h5 7_h5/OUTPUT/phelel_params.hdf5
```

当前 `8_cvq/extract_gcv.py` 默认读取：

```text
7_h5/OUTPUT/vaspelph.h5
```

如果这个文件不存在，说明第 7 步 `ELPH_DRIVER=mels` 还没有成功，不能继续算 `g_cvq`。

## 3. 提取 `g_cvq`

运行：

```bash
cd /Users/huhaitao/Desktop/code/Material/19.V2O5
/opt/miniconda3/bin/python 8_cvq/extract_gcv.py
```

它会生成：

```text
8_cvq/gcv_table.dat
8_cvq/gcv_summary.txt
```

`gcv_table.dat` 的列含义是：

```text
nu
phonon_meV
phonon_eV
Nph_300K
Re_g_eV
Im_g_eV
abs_g_eV
abs_g2_eV2
abs_weight_eV2 = |g|^2 * Nph
em_weight_eV2  = |g|^2 * (Nph + 1)
```

当前 V2O5 的关键检查结果应为：

```text
q [-0.5,0.5) = [-0.46511628 -0.5         0.        ]
q [0,1)      = [0.53488372 0.5        0.        ]
VBAND = 56
CBAND = 57
声子支数 = 42
```

如果 `q`、带号或声子支数不对，先不要继续算 `A_ind`。

## 4. 计算 `A_ind`

运行：

```bash
cd /Users/huhaitao/Desktop/code/Material/19.V2O5
/opt/miniconda3/bin/python 8_cvq/calc_A_ind.py
```

它会读取：

```text
8_cvq/gcv_table.dat
7_h5/POSCAR
5_meff/OUTPUT/meff    # 只用于可选读取 F0 候选值
```

并生成：

```text
8_cvq/A_ind_table.dat
8_cvq/A_ind_summary.txt
```

当前脚本采用的 `D_ph` 换算约定是：

```text
D_nu = |g_nu| * sqrt(2 * M_cell * epsilon_nu) / hbar
```

输出单位为：

```text
D_nu: eV/cm
epsilon_nu: eV
A_coeff_per_F0_5_2 = A_ind / F0^(5/2)
```

也就是说，在没有明确给定 `F0` 时，脚本给出的主结果是：

```text
A_ind = A_coeff_per_F0_5_2 * F0^(5/2)
```

当前 V2O5 的总结果是：

```text
A_ind,total = 6.919558033498e+15 * F0^(5/2)
```

如果临时把 `5_meff/OUTPUT/meff` 中的

```text
m_r_parallel/m0 = 1.160597
```

当作 `F0` 候选值代入，则：

```text
F0^(5/2) = 1.451120939869
A_ind,total = 1.004111555705e+16
```

但这里要注意：`F0` 的物理定义必须和你的模型公式一致。脚本不会自动判断 `F0` 是否应该等于沿场约化质量。

## 5. 为什么正式值要全支求和

当前最大 `|g|` 的声子支是：

```text
nu = 19
epsilon_ph = 66.714059766400 meV
|g| = 4.002518227319e-02 eV
A_nu19 = 1.734931980683e+15 * F0^(5/2)
占总和 25.07%
```

但对 `A_ind` 贡献最大的支是：

```text
nu = 5
epsilon_ph = 16.378678644590 meV
|g| = 2.819123061865e-02 eV
A_nu5 = 2.412031932964e+15 * F0^(5/2)
占总和 34.86%
```

原因是公式里有 `(1 + 2N_ph)`。低能声子的 Bose 占据数更大，所以最大 `|g|` 不一定等于最大 `A` 贡献。

因此建议：

```text
正式 A_ind: 用所有声子支求和
近似分析: 可以列出最大贡献支、最大 |g| 支、前几大贡献支
```

## 6. 换别的体系时需要改哪里

这两个脚本目前是为当前 V2O5 路径和跃迁写的，不是完全免修改的通用脚本。换体系时至少要检查并修改下面几处。

### 6.1 `extract_gcv.py`

文件：

```text
8_cvq/extract_gcv.py
```

需要改：

```python
H5_FILE = ROOT / "7_h5" / "OUTPUT" / "vaspelph.h5"
KV = np.array([0.4651162791, 0.5, 0.0])
KC = np.array([0.0, 0.0, 0.0])
VBAND = 56
CBAND = 57
SPIN = 0
TEMPERATURE = 300.0
```

含义：

```text
H5_FILE      你的 vaspelph.h5 路径
KV           初态价带 k 点，通常是 VBM
KC           末态导带 k 点，通常是 CBM
VBAND        初态价带的 VASP 1-based band index
CBAND        末态导带的 VASP 1-based band index
SPIN         非磁性通常为 0；自旋极化要确认取哪一个 spin channel
TEMPERATURE  Bose 因子温度，默认 300 K
```

换体系前必须重新确认 `VBAND/CBAND`。不要照抄 V2O5 的 56/57。

还要确认 `phonon_eigenvalues` 的单位。当前 VASP 输出在这个 h5 里是 meV 量级，所以脚本用了：

```python
phonon_mev = phonon_e[ikp_final, ik_initial, :]
eph = phonon_mev / 1000.0
```

如果你的 h5 中声子能量已经是 eV，需要把这里改成：

```python
eph = phonon_e[ikp_final, ik_initial, :]
phonon_mev = eph * 1000.0
```

判断方法：如果数值大约是 `10-200`，通常是 meV；如果大约是 `0.01-0.2`，通常是 eV。

### 6.2 `calc_A_ind.py`

文件：

```text
8_cvq/calc_A_ind.py
```

需要改：

```python
GC_TABLE = OUTDIR / "gcv_table.dat"
POSCAR = ROOT / "7_h5" / "POSCAR"
MEFF_REPORT = ROOT / "5_meff" / "OUTPUT" / "meff"
TEMPERATURE_K = 300.0
A_CONSTANT = 5.3464012e-4
ATOMIC_MASS_AMU = {
    "O": 15.9994,
    "V": 50.9415,
}
```

含义：

```text
GC_TABLE       上一步生成的 gcv_table.dat
POSCAR         对应 primitive cell 的 POSCAR
MEFF_REPORT    可选；如果你想让脚本读 m_r_parallel/m0 作为 F0 候选值
TEMPERATURE_K  要和 extract_gcv.py 中的温度一致
A_CONSTANT     你的公式前因子；如果模型公式变了，这里也要变
ATOMIC_MASS_AMU  体系元素的原子质量表
```

如果换成其他元素，必须在 `ATOMIC_MASS_AMU` 里补全所有元素。例如 MoS2：

```python
ATOMIC_MASS_AMU = {
    "Mo": 95.95,
    "S": 32.06,
}
```

当前脚本还有一个 V2O5 专用检查：

```python
if len(nu) != 42:
    raise ValueError(...)
```

换体系时要把 `42` 改成 `3 * 原胞原子数`，或者改成从 `POSCAR` 自动判断。比如原胞 3 个原子就是 `9` 支，原胞 20 个原子就是 `60` 支。

## 7. 这几个 py 文件能不能直接算别的体系

结论：不能无脑直接用，但可以作为模板快速改。

### 可以直接复用的部分

这些逻辑通常可以复用：

```text
读取 vaspelph.h5 的 kpoints/vkpt_k、kpoints/vkpt_kp
读取 matrix_elements/elph
读取 band_start_k、band_start_kp 并映射 VASP band index
按 kv、kc 找最近 k 点
提取所有 nu 的复数 g
计算 |g|、|g|^2、Nph
按所有声子支求和 A_ind
输出 table 和 summary
```

### 不能直接照搬的部分

这些必须按体系改：

```text
H5 文件路径
初态 k 点 KV
末态 k 点 KC
VBAND/CBAND
SPIN
原胞 POSCAR 路径
元素质量表
声子支数检查
F0 的来源和定义
声子能量单位判断
```

所以对另一个体系，推荐做法是：

```bash
cp -r 8_cvq /path/to/new_system/8_cvq
```

然后按上面的清单修改脚本参数，再运行：

```bash
cd /path/to/new_system
/opt/miniconda3/bin/python 8_cvq/extract_gcv.py
/opt/miniconda3/bin/python 8_cvq/calc_A_ind.py
```

## 8. 换体系的最小检查清单

1. `vaspelph.h5` 已经由正确的 `ELPH_DRIVER=mels` 计算生成。
2. `KPOINTS_ELPH` 里确实包含目标 `KV` 和 `KC`。
3. `VBAND` 和 `CBAND` 已从当前体系的 `OUTCAR/EIGENVAL` 确认。
4. `extract_gcv.py` 输出中的 `initial VBM k` 和 `final CBM k` 与目标一致。
5. `q` 与你要的间接跃迁一致。
6. 声子支数等于 `3N`。
7. `phonon_meV/phonon_eV` 单位判断正确。
8. `calc_A_ind.py` 的元素质量表覆盖 `POSCAR` 中所有元素。
9. `F0` 的定义已确认；如果没有确认，只报告 `A/F0^(5/2)`。
10. 最终 `A_ind_table.dat` 的 `frac_of_total` 求和应接近 1。

## 9. 当前结果文件

当前已经生成的结果在：

```text
8_cvq/gcv_table.dat
8_cvq/gcv_summary.txt
8_cvq/A_ind_table.dat
8_cvq/A_ind_summary.txt
```

查看摘要：

```bash
sed -n '1,40p' /Users/huhaitao/Desktop/code/Material/19.V2O5/8_cvq/gcv_summary.txt
sed -n '1,45p' /Users/huhaitao/Desktop/code/Material/19.V2O5/8_cvq/A_ind_summary.txt
```

重新计算：

```bash
cd /Users/huhaitao/Desktop/code/Material/19.V2O5
/opt/miniconda3/bin/python 8_cvq/extract_gcv.py
/opt/miniconda3/bin/python 8_cvq/calc_A_ind.py
```

