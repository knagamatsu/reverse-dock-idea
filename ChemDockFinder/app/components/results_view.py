import streamlit as st
import pandas as pd
from typing import Dict, Optional, List, Callable
import plotly.graph_objects as go
import plotly.express as px

def display_search_results(results_df: pd.DataFrame,
                         on_select: Optional[Callable] = None):
    """
    検索結果テーブルの表示
    
    Args:
        results_df: 結果のDataFrame
        on_select: 行選択時のコールバック関数
    """
    st.subheader("Search Results")
    
    if results_df.empty:
        st.info("No results found")
        return
    
    # 結果テーブルの表示
    selected_rows = st.data_editor(
        results_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "PDB ID": st.column_config.TextColumn(
                "PDB ID",
                help="Click to view structure",
                width="small",
            ),
            "Similarity": st.column_config.ProgressColumn(
                "Similarity (%)",
                format="%d",
                min_value=0,
                max_value=100,
            ),
            "Resolution": st.column_config.NumberColumn(
                "Resolution (Å)",
                format="%.1f",
            ),
        },
        key='protein_table'
    )
    
    # 行選択時のコールバック
    if on_select and selected_rows:
        on_select(selected_rows)

def display_docking_results(docking_results: Dict):
    """
    ドッキング結果の表示
    
    Args:
        docking_results: ドッキング計算の結果
    """
    st.subheader("Docking Results")
    
    if not docking_results:
        st.info("No docking results available")
        return
    
    # タブで結果を整理
    tabs = st.tabs(["Summary", "Energy Plot", "Interactions", "Details"])
    
    with tabs[0]:  # Summary
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Best Energy", f"{docking_results['best_energy']:.2f} kcal/mol")
        with col2:
            st.metric("Number of Poses", str(docking_results['num_poses']))
        
        # 結果の概要表
        st.markdown("### Result Summary")
        summary_df = pd.DataFrame({
            'Metric': ['Calculation ID', 'Status', 'Date', 'Total Runtime'],
            'Value': [
                docking_results.get('id', 'N/A'),
                docking_results.get('status', 'N/A'),
                docking_results.get('date', 'N/A'),
                docking_results.get('runtime', 'N/A')
            ]
        })
        st.table(summary_df)
    
    with tabs[1]:  # Energy Plot
        if 'energies' in docking_results:
            energies = docking_results['energies']
            
            # エネルギー分布のプロット
            fig = go.Figure()
            fig.add_trace(go.Box(
                y=energies,
                name='Binding Energy',
                boxpoints='all',
                jitter=0.3,
                pointpos=-1.8
            ))
            
            fig.update_layout(
                title='Distribution of Binding Energies',
                yaxis_title='Binding Energy (kcal/mol)',
                showlegend=False
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # エネルギーランキング
            st.markdown("### Energy Ranking")
            energy_df = pd.DataFrame({
                'Pose': range(1, len(energies) + 1),
                'Energy (kcal/mol)': energies
            })
            st.dataframe(energy_df)
    
    with tabs[2]:  # Interactions
        if 'interactions' in docking_results:
            interactions = docking_results['interactions']
            
            # 相互作用の表示
            st.markdown("### Key Interactions")
            int_df = pd.DataFrame(interactions)
            st.dataframe(int_df)
            
            # 相互作用の可視化
            if len(interactions) > 0:
                fig = px.bar(
                    int_df,
                    x='residue',
                    y='strength',
                    color='type',
                    title='Protein-Ligand Interactions'
                )
                st.plotly_chart(fig, use_container_width=True)
    
    with tabs[3]:  # Details
        st.markdown("### Calculation Parameters")
        # パラメータの表示
        if 'params' in docking_results:
            params = docking_results['params']
            for key, value in params.items():
                st.text(f"{key}: {value}")
        
        # 詳細な結果データ
        st.markdown("### Full Results")
        if 'full_results' in docking_results:
            st.json(docking_results['full_results'])

def create_download_buttons(docking_results: Dict):
    """
    結果のダウンロードボタンを作成
    
    Args:
        docking_results: ドッキング結果
    """
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("Download Results (CSV)"):
            # CSVファイルの作成とダウンロード
            results_df = pd.DataFrame({
                'Pose': range(1, len(docking_results['energies']) + 1),
                'Energy': docking_results['energies']
            })
            st.download_button(
                label="Download CSV",
                data=results_df.to_csv(index=False),
                file_name="docking_results.csv",
                mime="text/csv"
            )
    
    with col2:
        if st.button("Download Structure (PDB)"):
            # PDBファイルのダウンロード
            if 'best_pose' in docking_results:
                st.download_button(
                    label="Download PDB",
                    data=docking_results['best_pose'],
                    file_name="best_pose.pdb",
                    mime="chemical/x-pdb"
                )

def results_view_section(search_results: Optional[pd.DataFrame] = None,
                        docking_results: Optional[Dict] = None,
                        on_protein_select: Optional[Callable] = None):
    """
    結果表示セクション全体の表示
    
    Args:
        search_results: 検索結果のDataFrame
        docking_results: ドッキング結果
        on_protein_select: タンパク質選択時のコールバック
    """
    if search_results is not None:
        display_search_results(search_results, on_protein_select)
        
    if docking_results is not None:
        st.divider()
        display_docking_results(docking_results)
        create_download_buttons(docking_results)