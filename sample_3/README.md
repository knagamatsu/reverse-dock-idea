# Reverse Docking Tool

分子構造から相互作用する可能性のあるタンパク質を探索するための逆ドッキングツールです。直感的なWebインターフェースを通じて、分子設計からドッキングシミュレーションまでをシームレスに実行できます。

## 主な機能

- 分子構造エディタ (Ketcher) による化合物のデザイン
- SMILESコードによる分子入力
- 3D分子構造のリアルタイムプレビュー
- PDBデータベースからの類似構造検索
- タンパク質構造の3Dビューア
- カスタマイズ可能なドッキングシミュレーションパラメータ

## インストール方法

### 必要条件

- Python 3.8以上
- pip (Pythonパッケージマネージャー)

### セットアップ

1. リポジトリのクローン:
```bash
git clone <repository-url>
cd reverse-docking-tool
```

2. 必要なパッケージのインストール:
```bash
pip install -r requirements.txt
```

### 依存ライブラリ

- streamlit
- streamlit-ketcher
- py3Dmol
- rdkit
- pandas

## 使用方法

1. アプリケーションの起動:
```bash
streamlit run app.py
```

2. Webブラウザで `http://localhost:8501` にアクセス

### 基本的な操作手順

1. **分子構造の入力**
   - Ketcherエディタで直接描画
   - SMILESコードを入力
   - デフォルトのアミノメチルシクロヘキサンを使用

2. **類似構造の検索**
   - 入力された分子構造に基づいて自動的に検索
   - 類似度、解像度などで結果をフィルタリング

3. **ドッキングシミュレーション**
   - シミュレーションパラメータの設定
   - 「Run Docking」ボタンでシミュレーション開始

## シミュレーション設定

### 基本設定
- グリッドサイズ: 16-32Å (デフォルト: 24Å)
- 中心座標: X, Y, Z
- 探索の徹底度: 1-32 (デフォルト: 8)
- 結合モード数: 1-20 (デフォルト: 9)

### 高度な設定
- エネルギー範囲: 1-10 kcal/mol
- GPUアクセラレーション
- 水分子の考慮
- スコアリング関数の選択 (Vina, Vinardo, AD4)

## カスタマイズ

アプリケーションのスタイルはCSSを通じてカスタマイズ可能です。主な設定は以下の場所で行われています：

```python
st.markdown("""
    <style>
        .stApp > header {
            background-color: transparent;
        }
        # その他のスタイル設定
    </style>
""", unsafe_allow_html=True)
```

## データ構造

### PDBデータの形式
結果テーブルには以下の情報が含まれます：
- PDB ID
- 類似度 (%)
- 解像度 (Å)
- 生物種
- 実験手法

## 注意事項

- PDBファイルは `./data` ディレクトリに配置する必要があります
- 大きな分子構造の場合、3D構造の生成に時間がかかる場合があります
- GPUアクセラレーションを使用する場合は、適切なドライバとCUDAのインストールが必要です

## トラブルシューティング

1. 3D構造生成エラー
   - SMILES形式が正しいか確認
   - RDKitのバージョンを確認

2. PDB構造が表示されない
   - データディレクトリのパスを確認
   - ファイル形式が正しいか確認

## ライセンス

[ライセンス情報を記載]

## 貢献

プルリクエストやイシューの報告を歓迎します。大きな変更を加える場合は、まずイシューで変更内容を議論してください。