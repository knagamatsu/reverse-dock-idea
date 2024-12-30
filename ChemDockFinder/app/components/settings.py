import streamlit as st
from typing import Dict, Optional

DEFAULT_SETTINGS = {
    "grid_size": 24,
    "exhaustiveness": 8,
    "num_modes": 9,
    "energy_range": 3,
    "center_x": 0.0,
    "center_y": 0.0,
    "center_z": 0.0,
    "scoring": "Vina",
    "use_gpu": False,
    "consider_waters": False,
    "cpu_threads": 4
}

def initialize_settings():
    """設定の初期化"""
    if 'docking_settings' not in st.session_state:
        st.session_state.docking_settings = DEFAULT_SETTINGS.copy()

def save_settings(settings: Dict):
    """
    設定の保存
    
    Args:
        settings: 設定パラメータ
    """
    st.session_state.docking_settings = settings
    st.success("Settings saved successfully")

def show_settings_dialog():
    """設定ダイアログの表示"""
    st.write("Configure Docking Parameters")
    
    # タブで設定を整理
    basic_tab, advanced_tab, system_tab = st.tabs(["Basic", "Advanced", "System"])
    
    current_settings = st.session_state.docking_settings.copy()
    
    with basic_tab:
        st.subheader("Basic Settings")
        
        # グリッド設定
        st.markdown("#### Grid Settings")
        current_settings['grid_size'] = st.slider(
            "Grid Size (Å)",
            min_value=16,
            max_value=32,
            value=current_settings['grid_size'],
            help="Size of the search space"
        )
        
        # 探索設定
        st.markdown("#### Search Parameters")
        current_settings['exhaustiveness'] = st.slider(
            "Exhaustiveness",
            min_value=1,
            max_value=32,
            value=current_settings['exhaustiveness'],
            help="Higher values increase accuracy but take longer"
        )
        current_settings['num_modes'] = st.slider(
            "Number of Binding Modes",
            min_value=1,
            max_value=20,
            value=current_settings['num_modes'],
            help="Number of poses to generate"
        )
    
    with advanced_tab:
        st.subheader("Advanced Settings")
        
        # エネルギー設定
        st.markdown("#### Energy Parameters")
        current_settings['energy_range'] = st.slider(
            "Energy Range (kcal/mol)",
            min_value=1,
            max_value=10,
            value=current_settings['energy_range'],
            help="Maximum energy difference between best and worst pose"
        )
        
        # スコアリング関数
        current_settings['scoring'] = st.selectbox(
            "Scoring Function",
            ["Vina", "Vinardo", "AD4"],
            index=["Vina", "Vinardo", "AD4"].index(current_settings['scoring']),
            help="Choose scoring function for pose evaluation"
        )
        
        # 水分子の考慮
        current_settings['consider_waters'] = st.checkbox(
            "Consider Water Molecules",
            value=current_settings['consider_waters'],
            help="Include water molecules in calculation"
        )
    
    with system_tab:
        st.subheader("System Settings")
        
        # GPU設定
        current_settings['use_gpu'] = st.checkbox(
            "Use GPU Acceleration",
            value=current_settings['use_gpu'],
            help="Enable GPU acceleration if available"
        )
        
        # CPU設定
        current_settings['cpu_threads'] = st.number_input(
            "CPU Threads",
            min_value=1,
            max_value=32,
            value=current_settings['cpu_threads'],
            help="Number of CPU threads to use"
        )
    
    # 設定の保存
    if st.button("Apply Settings"):
        save_settings(current_settings)
        return current_settings
    
    return None

def display_current_settings():
    """現在の設定を表示"""
    settings = st.session_state.docking_settings
    
    st.markdown("#### Current Settings")
    st.markdown(f"""
    **Basic Parameters:**
    - Grid Size: {settings['grid_size']}Å
    - Exhaustiveness: {settings['exhaustiveness']}
    - Binding Modes: {settings['num_modes']}
    
    **Advanced Parameters:**
    - Energy Range: {settings['energy_range']} kcal/mol
    - Scoring: {settings['scoring']}
    - Consider Waters: {"Yes" if settings['consider_waters'] else "No"}
    
    **System Settings:**
    - GPU Acceleration: {"Enabled" if settings['use_gpu'] else "Disabled"}
    - CPU Threads: {settings['cpu_threads']}
    """)

def settings_section():
    """設定セクション全体の表示"""
    initialize_settings()
    
    st.subheader("Simulation Settings")
    
    # 設定表示/編集の選択
    action = st.radio(
        "Settings Action",
        ["View Current Settings", "Edit Settings"],
        horizontal=True
    )
    
    if action == "View Current Settings":
        display_current_settings()
    else:
        show_settings_dialog()

def get_current_settings() -> Dict:
    """
    現在の設定を取得
    
    Returns:
        Dict: 現在の設定パラメータ
    """
    initialize_settings()
    return st.session_state.docking_settings