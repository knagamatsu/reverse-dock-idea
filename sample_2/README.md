# AutoDock Vina 環境構築手順

以下の手順は、Ubuntu 24.04 LTS 上で AutoDock Vina を使用するための環境構築方法を示しています。

## 前提条件

* **OS**: Ubuntu 24.04 LTS
* **Python バージョン**: 3.9

## 環境構築手順

### 1. Anaconda のインストール

Anaconda をインストールしていない場合は、以下の手順でインストールしてください。

```bash
# Anaconda インストーラの実行
bash Anaconda3-*.sh

# conda コマンドの確認
conda --version
```

### 2. 仮想環境の作成と有効化

```bash
# 仮想環境の作成
conda create -n vina python=3.9

# 仮想環境の有効化
conda activate vina
```

### 3. 依存パッケージのインストール

```bash
# 必要なパッケージのインストール
conda install -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme
```

### 4. AutoDock Vina のインストール

```bash
# Vina のインストール
pip install vina
```

### 5. インストールの確認

```bash
# バージョン確認
python -c "import vina; print(vina.__version__)"

# ヘルプの表示
vina --help
```

### 6. サンプルデータのダウンロード

```bash
# 受容体とリガンドの構造データをダウンロード
wget https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/python_scripting/1iep_receptor.pdbqt
wget https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/example/python_scripting/1iep_ligand.pdbqt
```

### 7. ドッキングシミュレーションの実行

`vina_sample.py` を作成し、以下のコードを記述します：

```python
from vina import Vina

v = Vina(sf_name='vina')
v.set_receptor('1iep_receptor.pdbqt')
v.set_ligand_from_file('1iep_ligand.pdbqt')
v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
v.dock()
v.write_poses('1iep_ligand_out.pdbqt', n_poses=1)
```

実行結果の例：

```
Computing Vina grid ... done.
Performing docking (random seed: 753430492) ... 
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -13.26          0          0
   2        -11.3      3.011      12.42
   3       -11.15      3.801      12.29
   4       -10.18      2.555      12.57
   5       -9.774      2.944      12.56
   6       -8.984      3.931      12.67
   7       -8.805      3.531      12.12
   8       -7.426      7.414      12.61
   9       -7.407      7.359      11.96
  10       -7.403      7.358      11.99
  11       -7.381      7.155      12.54
  12       -7.075      7.443      12.35
  13       -6.904       4.18       12.6
  14       -6.564      1.551      2.175
  15       -6.393      7.175      12.46
  16       -6.088      5.357      12.87
  17       -5.628      3.283      12.41
  18       -5.569      3.751      13.21
  19       -5.485      3.889         13
  20       -2.511      4.055      14.15
```

## 参考文献

* [AutoDock Vina 1.2 のインストール - Hira Labo](リンク)
* [AutoDock Vina 1.2 のチュートリアル - Hira Labo](リンク)
* [AutoDock Vina 公式ドキュメント](リンク)