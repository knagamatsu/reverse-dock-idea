## 🛠 Development Setup

### Ubuntu環境でのセットアップ（推奨）
```bash
# 1. Minicondaのインストール
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

# 2. conda環境の作成
conda create -n vina_env python=3.10
conda activate vina_env
conda config --env --add channels conda-forge

# 3. 必要なパッケージのインストール
conda install -c conda-forge numpy swig boost-cpp libboost
conda install -c conda-forge rdkit

# 4. Python関連パッケージのインストール
pip install vina
pip install meeko
pip install streamlit streamlit-ketcher
pip install py3dmol

# 5. AutoDock Vinaバイナリのインストール
# GitHubからの最新リリースをダウンロード
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.3/vina_1.2.3_linux_x86_64
chmod +x vina_1.2.3_linux_x86_64
sudo mv vina_1.2.3_linux_x86_64 /usr/local/bin/vina
```

## 📊 Application Structure

```
chemistry-dashboard/
├── app.py                 # メインアプリケーション
├── requirements.txt       # 依存パッケージリスト
├── config/               # 設定ファイル
│   └── default.yaml     
├── modules/              # モジュール
│   ├── structure.py      # 構造処理モジュール
│   ├── docking.py       # ドッキング計算モジュール
│   └── visualization.py  # 可視化モジュール
├── data/                # データディレクトリ
│   ├── proteins/        # プリセットタンパク質
│   └── examples/        # サンプル分子
└── tests/              # テストコード
```

## 🔄 Development Log

### 2024-12-28
- プロジェクト開始
- 基本設計の策定
- 開発環境のセットアップ手順確立
- AutoDock Vinaをドッキングエンジンとして採用決定

### Key Decisions
1. **ドッキングエンジン選択**
   - AutoDock Vinaを採用
   - 理由：広いコミュニティ、充実したドキュメント、Pythonバインディング

2. **UI実装**
   - Streamlitを採用
   - 理由：高速な開発、化学構造エディタ(Ketcher)との統合が容易

3. **可視化**
   - py3Dmolを採用
   - 理由：軽量、インタラクティブ、Streamlitとの相性が良好

## 🎯 Next Steps
1. [ ] 基本的なStreamlitインターフェースの実装
2. [ ] Ketcherエディタの統合
3. [ ] AutoDock Vinaとの連携実装
4. [ ] 3D可視化機能の実装

## 🤝 Contribution Guidelines
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## 🐛 Known Issues
- Windows環境でのPythonバインディングの制限
- 大規模なタンパク質でのメモリ使用量の最適化が必要

## 🔗 References
- [AutoDock Vina Documentation](https://autodock-vina.readthedocs.io/)
- [Streamlit Documentation](https://docs.streamlit.io/)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [py3Dmol Documentation](https://3dmol.csb.pitt.edu/)