# ChemDockFinder

ChemDockFinderは、合成化合物の新規用途探索のための逆ドッキングツールです。化学者が手軽に利用できる、シンプルで実用的なインターフェースを提供します。

## 機能概要

- 分子構造エディタ（Ketcher）による直感的な構造入力
- SMILESコードによる分子指定
- 3D分子構造のリアルタイムプレビュー
- タンパク質構造の3D表示
- AutoDock Vinaによるブラインドドッキング
- 結果のインタラクティブな可視化

## インストール手順

### 1. 環境構築

```bash
# Condaによる環境作成
conda create -n vina_env python=3.9
conda activate vina_env
conda config --env --add channels conda-forge

# 基本パッケージのインストール
conda install -c conda-forge \
    rdkit \
    numpy \
    pandas \
    pathlib \
    swig \
    boost-cpp \
    openbabel

# Streamlit関連のインストール
pip install streamlit streamlit-ketcher

# 3D表示関連のインストール
pip install py3Dmol

# AutoDock Vina関連のインストール
pip install vina meeko
```

### 2. アプリケーションの取得

```bash
git clone https://github.com/yourusername/chemdockfinder.git
cd chemdockfinder
```

## 使用方法

```bash
# アプリケーションの起動
streamlit run app/main.py
```

ブラウザで自動的に http://localhost:8501 が開きます。

## インターフェース構成

アプリケーションは4つの主要なエリアで構成されています：

### 1. 構造入力エリア（左上）
- Ketcherエディタによる構造描画
- SMILES形式での入力
- デフォルト分子：アミノメチルシクロヘキサン（NCC1CCCCC1）

### 2. 検索結果エリア（右上）
- タンパク質のリスト表示
- PDB ID、類似度、解像度などの情報表示
- 実験手法の表示

### 3. 分子構造表示エリア（左下）
- 3D構造のリアルタイム表示
- ドッキングシミュレーション制御
- パラメータ設定機能

### 4. タンパク質構造表示エリア（右下）
- 選択したタンパク質の3D表示
- 結合部位の可視化

## プロジェクト構成

```
chemdockfinder/
├── app/
│   └── main.py           # Streamlitアプリケーション
├── core/
│   ├── docking.py        # Vinaドッキング処理
│   ├── structure.py      # 構造処理機能
│   └── database.py       # データベース管理
├── data/
│   └── proteins/         # タンパク質構造ファイル
└── requirements.txt
```

## シミュレーション設定

### 基本設定
- グリッドサイズ: 16-32Å (デフォルト: 24Å)
- 探索の徹底度: 1-32 (デフォルト: 8)
- 結合モード数: 1-20 (デフォルト: 9)
- エネルギー範囲: 1-10 kcal/mol

### 高度な設定
- GPUアクセラレーション
- 水分子の考慮
- スコアリング関数の選択（Vina, Vinardo, AD4）

## 開発状況

### 現在の実装状況
- [x] 基本的なUI実装
- [x] 分子構造入力（Ketcher）
- [x] 3D構造表示
- [x] シミュレーション設定ダイアログ
- [x] サンプルデータでの結果表示

### 開発予定
- [ ] AutoDock Vinaとの完全統合
- [ ] PDBデータベースとの連携
- [ ] 結果の保存機能
- [ ] バッチ処理対応

## 注意事項

- 大きな分子構造の場合、3D構造生成に時間がかかる場合があります
- ブラインドドッキングは計算時間が長くなる可能性があります
- タンパク質構造ファイルは事前に準備が必要です

## 技術仕様

- **フロントエンド**: Streamlit
- **分子表示**: py3Dmol
- **構造入力**: streamlit-ketcher
- **ドッキング**: AutoDock Vina
- **構造処理**: RDKit, OpenBabel

