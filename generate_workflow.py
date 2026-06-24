import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Define a color palette
COLOR_PRIMARY = '#1e3a8a'     # Deep Blue
COLOR_ACCENT = '#2563eb'      # Blue
COLOR_SUCCESS = '#16a34a'     # Green
COLOR_WARNING = '#ea580c'     # Orange/Yellow
COLOR_BG_COLUMN = '#f8fafc'   # Slate 50
COLOR_BORDER_COLUMN = '#cbd5e1' # Slate 300
COLOR_TEXT_MAIN = '#0f172a'   # Slate 900
COLOR_TEXT_MUTED = '#64748b'  # Slate 500

def create_workflow_diagram(output_path):
    fig, ax = plt.subplots(figsize=(16, 9.5), dpi=150)
    fig.patch.set_facecolor('#ffffff')
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # ----------------------------------------------------
    # Titles
    # ----------------------------------------------------
    ax.text(0.5, 0.96, "Joint miRNA-mRNA Network Analysis Pipeline", 
            fontsize=24, fontweight='bold', ha='center', va='center', color=COLOR_PRIMARY)
    ax.text(0.5, 0.92, "(Confounder Adjustment, Residual Correlation, Functional Annotation & ML Classification)", 
            fontsize=14, fontstyle='italic', ha='center', va='center', color=COLOR_TEXT_MUTED)

    # Helper function to draw columns
    def draw_column(x_start, width, title, num_items=None):
        rect = mpatches.FancyBboxPatch((x_start, 0.06), width, 0.82, boxstyle="round,pad=0.01,rounding_size=0.02",
                                      facecolor=COLOR_BG_COLUMN, edgecolor=COLOR_BORDER_COLUMN, linewidth=1.5, zorder=0)
        ax.add_patch(rect)
        ax.text(x_start + width/2, 0.85, title, fontsize=15, fontweight='bold', ha='center', va='center', color=COLOR_PRIMARY)

    # Draw columns
    draw_column(0.02, 0.16, "Input Data")
    draw_column(0.20, 0.60, "Workflow Phases")
    draw_column(0.82, 0.16, "Outputs")

    # ----------------------------------------------------
    # Helper to draw workflow steps/cards
    # ----------------------------------------------------
    def draw_card(title, x, y, w, h, step_num="", card_type="primary"):
        # Select color based on card type
        if card_type == "primary":
            border_c = '#3b82f6'
            bg_c = '#ffffff'
            num_c = '#2563eb'
        elif card_type == "ml":
            border_c = '#10b981'
            bg_c = '#ffffff'
            num_c = '#059669'
        elif card_type == "report":
            border_c = '#f97316'
            bg_c = '#fff7ed'
            num_c = '#ea580c'
        else:
            border_c = COLOR_BORDER_COLUMN
            bg_c = '#ffffff'
            num_c = COLOR_TEXT_MUTED

        # Draw card background
        card = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.005,rounding_size=0.015",
                                      facecolor=bg_c, edgecolor=border_c, linewidth=1.5, zorder=1)
        ax.add_patch(card)

        # Title
        ax.text(x + w/2, y + h - 0.02, title, fontsize=9.5, fontweight='bold', ha='center', va='center', color=COLOR_TEXT_MAIN, zorder=2)
        
        # Step Number Badge
        if step_num:
            ax.text(x + 0.015, y + h - 0.02, step_num, fontsize=10, fontweight='bold', ha='center', va='center', color=num_c, zorder=2)

    # ----------------------------------------------------
    # Column 1: Input Data Cards
    # ----------------------------------------------------
    # 1. Raw expression
    draw_card("Raw Transcript Data", 0.03, 0.64, 0.14, 0.16, card_type="input")
    ax.text(0.10, 0.73, "Raw Expression Matrices\n(miRNA & mRNA files)", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)
    # 2. Clinical Metadata
    draw_card("Clinical Metadata", 0.03, 0.40, 0.14, 0.16, card_type="input")
    ax.text(0.10, 0.49, "Donor Covariates\n(Diagnosis, RIN, PMI,\nAge, Sex)", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)
    # 3. Target database
    draw_card("Reference Databases", 0.03, 0.16, 0.14, 0.16, card_type="input")
    ax.text(0.10, 0.25, "Target Interactions\n(miRTarBase / TarBase\ncurated JSON mapping)", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)

    # ----------------------------------------------------
    # Column 2: Workflow Phase Cards
    # ----------------------------------------------------
    # Row 1: Steps 1, 2, 3
    # Step 1: Preprocessing & QC
    draw_card("Quality Control", 0.22, 0.62, 0.17, 0.18, step_num="01", card_type="primary")
    # Step 2: VIF Multicollinearity
    draw_card("Collinearity Check", 0.415, 0.62, 0.17, 0.18, step_num="02", card_type="primary")
    # Step 3: OLS Regress
    draw_card("Confounder Regress", 0.61, 0.62, 0.17, 0.18, step_num="03", card_type="primary")

    # Row 2: Steps 4, 5, 6
    # Step 4: Residual Pearson Correlation
    draw_card("Residual Correlation", 0.22, 0.38, 0.17, 0.18, step_num="04", card_type="primary")
    # Step 5: Target ORA
    draw_card("Target Enrichment", 0.415, 0.38, 0.17, 0.18, step_num="05", card_type="primary")
    # Step 6: ML Classification
    draw_card("Diagnostic ML", 0.61, 0.38, 0.17, 0.18, step_num="06", card_type="ml")

    # Row 3: Step 7 (Generate Report & Network Synthesis)
    draw_card("Integrative Report & Network Synthesis", 0.22, 0.12, 0.56, 0.18, step_num="07", card_type="report")

    # ----------------------------------------------------
    # Column 3: Outputs
    # ----------------------------------------------------
    draw_card("QC & EDA Visuals", 0.83, 0.70, 0.14, 0.11, card_type="output")
    ax.text(0.90, 0.73, "Preprocessed & Normalised\nBoxplots & Distributions", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)

    draw_card("OLS Coefficients", 0.83, 0.55, 0.14, 0.11, card_type="output")
    ax.text(0.90, 0.58, "Covariate-Adjusted Betas\n& Volcano Plots", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)

    draw_card("Robust Interactions", 0.83, 0.40, 0.14, 0.11, card_type="output")
    ax.text(0.90, 0.43, "miRNA-mRNA target pairs\n& Residual Scatter Plots", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)

    draw_card("Pathways & PPI", 0.83, 0.25, 0.14, 0.11, card_type="output")
    ax.text(0.90, 0.28, "Reactome Pathways\n& STRING PPI validation", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)

    draw_card("Classifier Metrics", 0.83, 0.10, 0.14, 0.11, card_type="output")
    ax.text(0.90, 0.13, "Diagnostic Performance\n(ROC-AUC scores)", fontsize=8, ha='center', va='center', color=COLOR_TEXT_MUTED)

    # ----------------------------------------------------
    # INSETS: Embedding mini-plots inside workflow cards!
    # ----------------------------------------------------
    # 01 QC: Mini Boxplot
    ax_ins1 = fig.add_axes([0.23, 0.63, 0.15, 0.11])
    ax_ins1.axis('off')
    # Draw simple boxes
    ax_ins1.add_patch(mpatches.Rectangle((1, 0.2), 0.8, 0.6, facecolor='#93c5fd', edgecolor='#2563eb', alpha=0.7))
    ax_ins1.add_patch(mpatches.Rectangle((3, 0.3), 0.8, 0.5, facecolor='#cbd5e1', edgecolor='#475569', alpha=0.7))
    ax_ins1.plot([1.4, 1.4], [0.1, 0.9], color='#2563eb', linewidth=1.5)
    ax_ins1.plot([3.4, 3.4], [0.1, 0.9], color='#475569', linewidth=1.5)
    ax_ins1.plot([1.1, 1.7], [0.5, 0.5], color='#2563eb', linewidth=2)
    ax_ins1.plot([3.1, 3.7], [0.55, 0.55], color='#475569', linewidth=2)
    ax_ins1.set_xlim(0, 5)
    ax_ins1.set_ylim(0, 1)

    # 02 VIF: Mini Bar Chart
    ax_ins2 = fig.add_axes([0.425, 0.63, 0.15, 0.11])
    ax_ins2.axis('off')
    bars = ax_ins2.bar(['A', 'P', 'R', 'S'], [1.85, 1.30, 1.43, 1.82], color='#3b82f6', alpha=0.8, width=0.6)
    ax_ins2.axhline(y=5.0, color='#ef4444', linestyle='--', linewidth=1, alpha=0.7)
    ax_ins2.set_ylim(0, 6)
    # Add tiny text labels
    for bar in bars:
        height = bar.get_height()
        ax_ins2.text(bar.get_x() + bar.get_width()/2., height + 0.2, f'{height:.1f}',
                    ha='center', va='bottom', fontsize=7, color=COLOR_TEXT_MUTED)

    # 03 OLS: Mini formula / regression plot
    ax_ins3 = fig.add_axes([0.62, 0.63, 0.15, 0.11])
    ax_ins3.axis('off')
    ax_ins3.text(0.5, 0.7, r"$\log_2(Y) = X\beta + \varepsilon$", fontsize=10, ha='center', va='center', color=COLOR_PRIMARY)
    # Draw simple regression line
    xs = np.linspace(0.1, 0.9, 10)
    ys = 0.8 - 0.6 * xs + np.random.normal(0, 0.05, 10)
    ax_ins3.scatter(xs, ys, color='#64748b', s=8, alpha=0.6)
    ax_ins3.plot(xs, 0.8 - 0.6 * xs, color='#ef4444', linewidth=1.5)
    ax_ins3.set_xlim(0, 1)
    ax_ins3.set_ylim(0, 1)

    # 04 Residual Correlation: Scatter Plot
    ax_ins4 = fig.add_axes([0.23, 0.39, 0.15, 0.11])
    ax_ins4.axis('off')
    # Generate anti-correlated scatter
    np.random.seed(42)
    rx = np.random.uniform(0.1, 0.9, 20)
    ry = 0.9 - 0.8 * rx + np.random.normal(0, 0.08, 20)
    ax_ins4.scatter(rx, ry, color='#3b82f6', s=8, alpha=0.7)
    ax_ins4.plot(np.unique(rx), np.poly1d(np.polyfit(rx, ry, 1))(np.unique(rx)), color='#2563eb', linewidth=1.5)
    ax_ins4.text(0.8, 0.8, "r = -0.91", fontsize=8.5, weight='bold', color='#1e3a8a', ha='right')
    ax_ins4.set_xlim(0, 1)
    ax_ins4.set_ylim(0, 1)

    # 05 Target ORA: Enrichment bars
    ax_ins5 = fig.add_axes([0.425, 0.39, 0.15, 0.11])
    ax_ins5.axis('off')
    pathways = ['miR-9', 'miR-132', 'miR-155']
    p_vals = [2.7, 1.8, 1.2]  # -log10 p-values
    y_pos = np.arange(len(pathways))
    ax_ins5.barh(y_pos, p_vals, color='#10b981', alpha=0.8, height=0.5)
    ax_ins5.axvline(x=1.3, color='#ef4444', linestyle='--', linewidth=0.8, alpha=0.7) # p=0.05
    ax_ins5.set_yticks(y_pos)
    ax_ins5.set_yticklabels(pathways, fontsize=6.5, color=COLOR_TEXT_MAIN)
    ax_ins5.tick_params(left=False, bottom=False, labelbottom=False)

    # 06 Diagnostic ML: ROC Curve
    ax_ins6 = fig.add_axes([0.62, 0.39, 0.15, 0.11])
    ax_ins6.axis('off')
    # ROC curve
    roc_x = np.linspace(0, 1, 100)
    roc_y = 1 - (1 - roc_x)**3 # high AUC
    ax_ins6.plot(roc_x, roc_y, color='#10b981', linewidth=2)
    ax_ins6.plot([0, 1], [0, 1], color='#cbd5e1', linestyle='--', linewidth=1)
    ax_ins6.text(0.65, 0.25, "AUC = 0.938", fontsize=8.5, weight='bold', color='#047857', ha='center')
    ax_ins6.set_xlim(0, 1)
    ax_ins6.set_ylim(0, 1)

    # 07 Report & Network: Complex network plot + report icons
    ax_ins7 = fig.add_axes([0.23, 0.13, 0.54, 0.11])
    ax_ins7.axis('off')
    # Draw simple report items
    # Doc icon 1: Summary Stats
    ax_ins7.add_patch(mpatches.Rectangle((0.05, 0.15), 0.08, 0.7, facecolor='#ffffff', edgecolor='#ea580c', linewidth=1))
    ax_ins7.text(0.09, 0.5, "Summary\nStats", fontsize=8, ha='center', va='center')
    # Doc icon 2: Plots
    ax_ins7.add_patch(mpatches.Rectangle((0.20, 0.15), 0.08, 0.7, facecolor='#ffffff', edgecolor='#ea580c', linewidth=1))
    ax_ins7.text(0.24, 0.5, "Volcano\nPlots", fontsize=8, ha='center', va='center')
    # Doc icon 3: Pathways
    ax_ins7.add_patch(mpatches.Rectangle((0.35, 0.15), 0.08, 0.7, facecolor='#ffffff', edgecolor='#ea580c', linewidth=1))
    ax_ins7.text(0.39, 0.5, "Reactome\nPaths", fontsize=8, ha='center', va='center')
    
    # Draw simple Network Representation in the remaining space of Step 7
    # Nodes:
    nodes = [(0.65, 0.5), (0.75, 0.7), (0.85, 0.6), (0.72, 0.3), (0.82, 0.4), (0.92, 0.5)]
    labels = ['miR-132', 'EP300', 'FOXO1', 'SIRT1', 'miR-9', 'MAPK1']
    node_colors = ['#f43f5e', '#60a5fa', '#34d399', '#34d399', '#f43f5e', '#60a5fa']
    
    # Edges:
    edges = [(0, 1), (0, 2), (4, 3), (3, 2), (5, 2), (1, 2)]
    for e in edges:
        n1, n2 = nodes[e[0]], nodes[e[1]]
        ax_ins7.plot([n1[0], n2[0]], [n1[1], n2[1]], color='#cbd5e1', linewidth=1, zorder=1)
        
    for i, (nx, ny) in enumerate(nodes):
        ax_ins7.scatter(nx, ny, color=node_colors[i], s=70, edgecolor='#475569', zorder=2)
        ax_ins7.text(nx, ny - 0.14, labels[i], fontsize=7, weight='bold', ha='center', color=COLOR_TEXT_MAIN, zorder=3)
    ax_ins7.set_xlim(0, 1)
    ax_ins7.set_ylim(0, 1)

    # ----------------------------------------------------
    # CONNECTING ARROWS & FLOWS
    # ----------------------------------------------------
    # Input Data -> QC Card
    ax.annotate('', xy=(0.21, 0.71), xytext=(0.175, 0.71),
                arrowprops=dict(arrowstyle="->", color=COLOR_ACCENT, lw=2.5))
    
    # Input Data -> Clinical Metadata Metadata -> QC Card
    # Handled by metadata inputs going into step 1 and 2
    # Arrow between Row 1 Cards: Step 01 -> Step 02
    ax.annotate('', xy=(0.405, 0.71), xytext=(0.395, 0.71),
                arrowprops=dict(arrowstyle="->", color=COLOR_ACCENT, lw=2.5))
    # Step 02 -> Step 03
    ax.annotate('', xy=(0.60, 0.71), xytext=(0.59, 0.71),
                arrowprops=dict(arrowstyle="->", color=COLOR_ACCENT, lw=2.5))
    
    # Step 03 -> Step 04 (Adjusted expression goes to correlation)
    # Down arrow and left arrow
    ax.annotate('', xy=(0.305, 0.57), xytext=(0.695, 0.61),
                arrowprops=dict(arrowstyle="->", color=COLOR_ACCENT, lw=2.0, connectionstyle="arc3,rad=0.15"))
    
    # Step 03 -> Step 06 (ML classification on adjusted residuals)
    # Down arrow
    ax.annotate('', xy=(0.695, 0.57), xytext=(0.695, 0.61),
                arrowprops=dict(arrowstyle="->", color=COLOR_SUCCESS, lw=2.0))

    # Step 04 -> Step 05 (Residuals to target ORA)
    ax.annotate('', xy=(0.405, 0.47), xytext=(0.395, 0.47),
                arrowprops=dict(arrowstyle="->", color=COLOR_ACCENT, lw=2.5))

    # Step 05 -> Step 07 (Enrichments to Integrative report)
    ax.annotate('', xy=(0.49, 0.31), xytext=(0.49, 0.37),
                arrowprops=dict(arrowstyle="<-", color=COLOR_WARNING, lw=2.0, linestyle='--'))

    # Step 06 -> Step 07 (ML classification to Integrative report)
    ax.annotate('', xy=(0.695, 0.31), xytext=(0.695, 0.37),
                arrowprops=dict(arrowstyle="<-", color=COLOR_WARNING, lw=2.0, linestyle='--'))

    # Step 04 -> Step 07 (Residual correlation to Integrative report)
    ax.annotate('', xy=(0.305, 0.31), xytext=(0.305, 0.37),
                arrowprops=dict(arrowstyle="<-", color=COLOR_WARNING, lw=2.0, linestyle='--'))

    # ----------------------------------------------------
    # Workflow Phases to Outputs Column Arrows
    # ----------------------------------------------------
    # Connect Step 01 -> Output 1
    ax.annotate('', xy=(0.825, 0.75), xytext=(0.395, 0.79),
                arrowprops=dict(arrowstyle="->", color=COLOR_TEXT_MUTED, lw=1.2, ls=':'))
    # Connect Step 03 -> Output 2
    ax.annotate('', xy=(0.825, 0.60), xytext=(0.785, 0.70),
                arrowprops=dict(arrowstyle="->", color=COLOR_TEXT_MUTED, lw=1.2, ls=':'))
    # Connect Step 04 -> Output 3
    ax.annotate('', xy=(0.825, 0.45), xytext=(0.395, 0.37),
                arrowprops=dict(arrowstyle="->", color=COLOR_TEXT_MUTED, lw=1.2, ls=':'))
    # Connect Step 07 -> Output 4
    ax.annotate('', xy=(0.825, 0.30), xytext=(0.785, 0.20),
                arrowprops=dict(arrowstyle="->", color=COLOR_WARNING, lw=1.2, ls='--'))
    # Connect Step 06 -> Output 5
    ax.annotate('', xy=(0.825, 0.15), xytext=(0.785, 0.37),
                arrowprops=dict(arrowstyle="->", color=COLOR_SUCCESS, lw=1.2, ls='--'))

    # ----------------------------------------------------
    # Legend at the Bottom
    # ----------------------------------------------------
    legend_y = 0.02
    ax.text(0.20, legend_y, "Legend:", fontsize=10, fontweight='bold', color=COLOR_TEXT_MAIN, ha='left', va='center')
    
    # Legend Items
    # 1. Primary Workflow
    ax.plot([0.28, 0.32], [legend_y, legend_y], color=COLOR_ACCENT, lw=2.5)
    ax.text(0.33, legend_y, "Primary Workflow", fontsize=9.5, color=COLOR_TEXT_MAIN, ha='left', va='center')
    
    # 2. Downstream Analysis (ML)
    ax.plot([0.48, 0.52], [legend_y, legend_y], color=COLOR_SUCCESS, lw=2.0)
    ax.text(0.53, legend_y, "Downstream ML", fontsize=9.5, color=COLOR_TEXT_MAIN, ha='left', va='center')
    
    # 3. Reporting / Validation Flow
    ax.plot([0.65, 0.69], [legend_y, legend_y], color=COLOR_WARNING, lw=2.0, linestyle='--')
    ax.text(0.70, legend_y, "Reporting Flow", fontsize=9.5, color=COLOR_TEXT_MAIN, ha='left', va='center')

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight', pad_inches=0.1)
    plt.close()
    print(f"Successfully generated workflow diagram at {output_path}")

if __name__ == "__main__":
    import sys
    path = "docs/figures/workflow.png"
    if len(sys.argv) > 1:
        path = sys.argv[1]
    create_workflow_diagram(path)
