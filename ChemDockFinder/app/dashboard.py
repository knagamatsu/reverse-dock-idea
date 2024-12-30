import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, timedelta
from typing import Dict, List, Optional
import pathlib

def display_summary_metrics(results_db):
    """サマリーメトリクスの表示"""
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        total_calcs = len(results_db.get_all_calculations())
        st.metric("Total Calculations", total_calcs)

    with col2:
        recent_calcs = len(results_db.get_recent_calculations(days=1))
        st.metric("Last 24 Hours", recent_calcs)

    with col3:
        success_rate = 95  # 実際のデータから計算する必要あり
        st.metric("Success Rate", f"{success_rate}%")

    with col4:
        avg_time = "2.5 min"  # 実際のデータから計算する必要あり
        st.metric("Avg. Calculation Time", avg_time)

def plot_calculation_history(results_db):
    """計算履歴のプロット"""
    df = results_db.get_all_calculations()
    if not df.empty:
        fig = px.bar(
            df,
            x='Date',
            title='Calculation History',
            height=300
        )
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("No calculation history available")

def display_recent_results(results_db):
    """最近の計算結果の表示"""
    st.subheader("Recent Calculations")
    
    recent_df = results_db.get_recent_calculations(limit=5)
    if not recent_df.empty:
        st.dataframe(
            recent_df,
            use_container_width=True,
            column_config={
                "Date": st.column_config.DatetimeColumn(
                    "Date",
                    format="DD/MM/YY HH:mm"
                ),
                "Status": st.column_config.TextColumn(
                    "Status",
                    help="Calculation status"
                )
            }
        )
    else:
        st.info("No recent calculations")

def plot_energy_distribution(results_db):
    """エネルギー分布のプロット"""
    st.subheader("Energy Distribution")
    
    # サンプルデータ（実際のデータベースから取得する必要あり）
    energies = {
        'range': ['< -10', '-10 to -8', '-8 to -6', '-6 to -4', '> -4'],
        'count': [10, 25, 35, 20, 10]
    }
    
    fig = go.Figure(data=[
        go.Bar(
            x=energies['range'],
            y=energies['count'],
            marker_color='rgb(55, 83, 109)'
        )
    ])
    
    fig.update_layout(
        title='Distribution of Binding Energies',
        xaxis_title='Energy Range (kcal/mol)',
        yaxis_title='Count',
        height=300
    )
    
    st.plotly_chart(fig, use_container_width=True)

def show_resource_usage():
    """システムリソース使用状況の表示"""
    st.subheader("System Resources")
    
    # サンプルデータ（実際のシステムから取得する必要あり）
    col1, col2 = st.columns(2)
    
    with col1:
        cpu_usage = 45  # %
        st.progress(cpu_usage / 100, "CPU Usage: {}%".format(cpu_usage))
        
    with col2:
        memory_usage = 60  # %
        st.progress(memory_usage / 100, "Memory Usage: {}%".format(memory_usage))

def show_active_calculations():
    """アクティブな計算の表示"""
    st.subheader("Active Calculations")
    
    # サンプルデータ（実際のデータベースから取得する必要あり）
    active_calcs = [
        {"id": "calc_001", "progress": 75, "eta": "2 min"},
        {"id": "calc_002", "progress": 30, "eta": "5 min"}
    ]
    
    if active_calcs:
        for calc in active_calcs:
            st.markdown(f"**Calculation ID:** {calc['id']}")
            st.progress(calc['progress'] / 100, f"Progress: {calc['progress']}% (ETA: {calc['eta']})")
    else:
        st.info("No active calculations")

def dashboard_page(results_db):
    """ダッシュボードページ全体の表示"""
    st.title("ChemDockFinder Dashboard")
    
    # サマリーメトリクス
    st.markdown("### Overview")
    display_summary_metrics(results_db)
    
    # 2列レイアウト
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # 計算履歴
        plot_calculation_history(results_db)
        
        # エネルギー分布
        plot_energy_distribution(results_db)
    
    with col2:
        # アクティブな計算
        show_active_calculations()
        
        # システムリソース
        show_resource_usage()
    
    # 最近の結果
    display_recent_results(results_db)

def mini_dashboard():
    """ミニダッシュボードの表示（サイドバー用）"""
    st.markdown("### Quick Stats")
    
    # サンプルデータ
    stats = {
        "Active": 2,
        "Completed Today": 15,
        "Queue": 3
    }
    
    for label, value in stats.items():
        st.metric(label, value)
    
    # 簡易システム状態
    st.markdown("### System Status")
    status = "🟢 Operational"  # または "🔴 Issues Detected"
    st.markdown(f"**Status:** {status}")

class DashboardManager:
    """ダッシュボード管理クラス"""
    
    def __init__(self, results_db, update_interval: int = 5):
        """
        初期化
        
        Args:
            results_db: 結果データベース
            update_interval: 更新間隔（秒）
        """
        self.results_db = results_db
        self.update_interval = update_interval
        self.last_update = None
    
    def initialize_session_state(self):
        """セッション状態の初期化"""
        if 'dashboard_data' not in st.session_state:
            st.session_state.dashboard_data = {
                'active_calculations': [],
                'system_resources': {
                    'cpu_usage': 0,
                    'memory_usage': 0
                },
                'last_update': None
            }
    
    def update_dashboard_data(self):
        """ダッシュボードデータの更新"""
        current_time = datetime.now()
        
        if (self.last_update is None or 
            (current_time - self.last_update).seconds >= self.update_interval):
            
            # データの更新（実際のデータ取得処理に置き換える）
            st.session_state.dashboard_data = {
                'active_calculations': self.results_db.get_active_calculations(),
                'system_resources': self._get_system_resources(),
                'last_update': current_time
            }
            
            self.last_update = current_time
    
    def _get_system_resources(self) -> Dict:
        """システムリソース情報の取得"""
        # 実際のシステムモニタリング処理に置き換える
        return {
            'cpu_usage': 45,
            'memory_usage': 60
        }
    
    def display_dashboard(self):
        """ダッシュボードの表示"""
        self.initialize_session_state()
        self.update_dashboard_data()
        dashboard_page(self.results_db)