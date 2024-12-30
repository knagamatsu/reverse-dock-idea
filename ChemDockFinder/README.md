# ChemDockFinder (Sample 4)

実用的な逆ドッキングツールの実装バージョンです。ダッシュボード機能を備え、AutoDock Vinaによる実際のドッキング計算と結果の可視化を行います。

## ディレクトリ構成

```
sample_4/
├── app/
│   ├── main.py                 # メインのStreamlitアプリ
│   ├── dashboard.py            # ダッシュボード機能
│   └── components/
│       ├── structure_input.py  # 構造入力関連
│       ├── results_view.py     # 結果表示関連
│       ├── ligand_view.py      # リガンド3D表示
│       ├── protein_view.py     # タンパク質構造表示
│       └── settings.py         # 設定ダイアログ
├── core/
│   ├── docking/
│   │   ├── vina.py            # Vina関連の処理
│   │   ├── converter.py       # 構造変換
│   │   └── calculator.py      # ドッキング計算
│   ├── structure/
│   │   ├── processor.py       # 構造処理
│   │   └── validator.py       # 入力検証
│   └── database/
│       ├── protein_db.py      # タンパク質DB管理
│       └── results_db.py      # 結果保存管理
├── data/
│   ├── proteins/              # タンパク質構造ファイル
│   │   ├── prepared/         # 準備済みの構造
│   │   └── original/         # オリジナルのPDB
│   ├── results/              # 計算結果
│   │   ├── docking/         # ドッキング結果
│   │   └── logs/            # 計算ログ
│   └── configs/              # 設定ファイル
├── tests/                    # テストコード
├── notebooks/                # 開発用ノートブック
├── requirements.txt          # 依存パッケージ
├── setup.py                  # セットアップスクリプト
└── README.md                 # このファイル
```

## 主な機能

### 1. インターフェース
- 4エリア構成のメイン画面
  - 構造入力（Ketcher）
  - 検索結果表示
  - リガンド3D表示
  - ドッキング結果の3D表示
- ダッシュボード機能
  - 計算状況の表示
  - 結果の統計情報
  - バッチ処理の管理

### 2. ドッキング機能
- AutoDock Vinaによる計算
- ブラインドドッキング対応
- バッチ処理サポート
- 結果の自動保存

### 3. 結果表示
- ドッキング後の構造可視化
- 結合エネルギーの表示
- 相互作用の分析
- 結果のエクスポート

## セットアップ

```bash
# 1. 環境作成
conda create -n vina_env python=3.9
conda activate vina_env
conda config --env --add channels conda-forge

# 2. 基本パッケージ
conda install -c conda-forge \
    rdkit \
    numpy \
    pandas \
    pathlib \
    swig \
    boost-cpp \
    openbabel

# 3. Streamlit関連
pip install streamlit streamlit-ketcher

# 4. 3D表示関連
pip install py3Dmol

# 5. AutoDock Vina関連
pip install vina meeko
```

## 使用方法

```bash
# アプリケーションの起動
streamlit run app/main.py
```

## 開発ステータス

### 実装済み機能
- [ ] 基本的なUI構成
- [ ] Vina計算エンジンの統合
- [ ] 構造入力機能
- [ ] 3D表示機能
- [ ] ダッシュボード基本機能

### 開発予定
- [ ] バッチ処理機能
- [ ] 結果の分析機能
- [ ] エクスポート機能
- [ ] プログレス表示の改善

## 注意事項

- タンパク質構造は事前に準備が必要です
- 大規模な計算にはリソースの考慮が必要です
- 結果は自動的に保存されます

## 参考文献

- AutoDock Vina Documentation
- Streamlit Documentation
- RDKit Documentation
- py3Dmol Documentation

## ライセンス

[ライセンス情報を記載]