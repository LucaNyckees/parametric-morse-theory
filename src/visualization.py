import plotly.graph_objects as go
from plotly.subplots import make_subplots


def plotly_persistence_diagram(coordinates: list) -> None:

    maximum = max([trio[2] for trio in coordinates])

    fig = make_subplots(rows=2, cols=2, subplot_titles=("0-cells", "1-cells", "2-cells", "3-cells"),
                        shared_xaxes=True, shared_yaxes=True)

    for dim in range(4):
        coord_dim = [trio for trio in coordinates if len(trio[0]) == dim + 1]
        births_dim = [trio[1] for trio in coord_dim]
        deaths_dim = [trio[2] for trio in coord_dim]
        z = list(range(maximum + 1))

        scatter_trace = go.Scatter(x=births_dim, y=deaths_dim, mode='markers',
                                   marker=dict(color="royalblue"), name=f"{dim}-cells")
        fig.add_trace(scatter_trace, row=(dim // 2) + 1, col=(dim % 2) + 1)

        diagonal_trace = go.Scatter(x=z, y=z, mode='lines', line=dict(color="dodgerblue"), showlegend=False)
        fig.add_trace(diagonal_trace, row=(dim // 2) + 1, col=(dim % 2) + 1)

        for i in range(len(coord_dim)):
            line_trace = go.Scatter(x=[births_dim[i], births_dim[i]], y=[births_dim[i], deaths_dim[i]],
                                    mode='lines', line=dict(color="darkviolet", dash="dot"), showlegend=False)
            fig.add_trace(line_trace, row=(dim // 2) + 1, col=(dim % 2) + 1)

    fig.update_layout(title_text="Life coordinates of critical cells", height=600, width=800)
    fig.update_xaxes(title_text="Births", row=2, col=1)
    fig.update_yaxes(title_text="Deaths", row=2, col=1)
    fig.update_yaxes(title_text="Deaths", row=1, col=1)
    fig.update_xaxes(title_text="Births", row=2, col=2)
    fig.show()
