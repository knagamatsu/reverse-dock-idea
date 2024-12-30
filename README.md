# Reverse Docking Tool Development

このリポジトリは分子構造から相互作用する可能性のあるタンパク質を探索するための逆ドッキングツールの開発プロセスを記録したものです。

## 📁 Repository Structure

```
reverse-docking-tool/
├── sample_1/              # 初期実装
│   ├── app.py            # Streamlitベースの実装
│   └── README.md         # 詳細な実装ドキュメント
├── sample_3/              # AutoDock Vina統合版
│   ├── app.py            # 新UI実装
│   └── README.md         # セットアップと利用ガイド
└── README.md             # このファイル
```

## 🔄 Development History

### Sample 1: 初期実装
- Streamlitベースの基本的なUI実装
- Ketcherエディタによる分子構造入力
- py3Dmolによる3D構造可視化
- 基本的なPDB検索機能

実装の詳細は [sample_1/README.md](./sample_1/README.md) を参照してください。

### Sample 3: AutoDock Vina統合版
- conda環境を使用した開発環境の整備
- AutoDock Vinaをドッキングエンジンとして統合
- モジュール構造の改善
- 高度なシミュレーション設定の実装

セットアップと使用方法の詳細は [sample_3/README.md](./sample_3/README.md) を参照してください。

## 💡 Implementation Comparison

### Sample 1
- **長所**:
  - シンプルな実装
  - 最小限の依存関係
  - 素早いプロトタイピング
- **短所**:
  - 限定的な機能
  - スケーラビリティの課題
  - ドッキング計算機能なし

### Sample 3
- **長所**:
  - 本格的なドッキング計算機能
  - モジュール化された設計
  - 拡張性の高い構造
- **短所**:
  - 複雑な環境設定
  - より多くの依存関係
  - 高いシステム要件

## 🎯 Future Development

### 短期目標
1. [ ] パフォーマンス最適化
2. [ ] エラーハンドリングの改善
3. [ ] ユーザーインターフェースの統一

### 長期目標
1. [ ] 機械学習モデルの統合
2. [ ] バッチ処理機能の実装
3. [ ] クラウドデプロイメント対応

## 🛠 Overall Technology Stack

- **フロントエンド**: Streamlit, Ketcher
- **バックエンド**: Python, AutoDock Vina
- **可視化**: py3Dmol, RDKit
- **開発環境**: conda, pip

## 📝 Notes

- 各サンプルディレクトリには独自のREADMEがあり、詳細な実装情報が記載されています
- 開発の進行に伴い、新しいアプローチやテクノロジーを試験的に導入しています
- 各実装は独立して動作可能です

## 🤝 Contributing

1. 適切なサンプルディレクトリを選択
2. 機能改善やバグ修正を実装
3. テストを追加
4. プルリクエストを作成

## 📚 References

- [Streamlit Documentation](https://docs.streamlit.io/)
- [AutoDock Vina Documentation](https://autodock-vina.readthedocs.io/)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [py3Dmol Documentation](https://3dmol.csb.pitt.edu/)