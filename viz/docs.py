import dash
from dash import html
import dash_bootstrap_components as dbc

from viz.web import make_tooltip, make_data_redirect_buttons, wrap_icon


def home_welcome_info():
    return [
        html.Img(src=dash.get_asset_url("website_header.png"), style={"width": "55%"}, alt="ContactTracing Header"),
        html.Br(),
        html.P(),
        html.Div(["Chromosomal Instability (CIN) is a hallmark of human cancer that is associated with metastasis and "
                  "immune evasion; yet how CIN shapes the tumor ecosystem remained poorly understood. Exploiting "
                  "intrinsic biological variability in single cell data, we developed a systems level approach called ",
                  html.I("ContactTracing"),
                  " that predicts ",
                  html.I(html.U("the effect")),
                  " of ligand-receptor-mediated interactions on the tumor microenvironment. Coupling ",
                  html.I("ContactTracing"),
                  " with isogenic breast cancer models distinguished by tumor cell-intrinsic rates of chromosomal "
                  "instability (CIN) or downstream activation of the cGAS-STING pathway, we identified tumor ligands "
                  "emanating from an ER-stress response as critical mediators of immune suppression in tumors with CIN. "
                  "Use this dashboard to explore how CIN-induced STING signaling in cancer cells shapes the tumor "
                  "microenvironment (TME). Or use ",
                  html.I("ContactTracing"),
                  " on your own data by checking out our ",
                  html.A("GitHub repository", href="https://github.com/LaughneyLab/ContactTracing_tutorial"),
                  "!",
                  html.P(),
                  html.Div([
                      dbc.Button(wrap_icon("fa-brands fa-github", "Run ", html.I("ContactTracing"), " on your data"),
                                 color="info", size='lg',
                                 href='https://github.com/LaughneyLab/ContactTracing_tutorial')
                  ], className='text-center')
              ]),
        html.P(),
        html.Div([
            html.H5("Key Features of ContactTracing"),
            html.Ul([
                html.Li("Identifies conditionally-dependent cell-cell interactions in the microenvironment."),
                html.Li("Infers ligand effects on gene expression in target (receptor-expressing) cells without prior "
                        "knowledge of downstream signaling."),
                html.Li("Captures cell type-specific responses to ligand-receptor-mediated interactions."),
                html.Li("Prioritizes responsive cell-cell interactions that more frequently co-localize in matched "
                        "spatial transcriptomics data.")
            ])
        ]),
        html.P(),
        html.Div([
            html.H5("Included Visualizations"),
            make_data_redirect_buttons()
        ]),
        html.P(),
        html.Div([html.H5("Citation"),
                dbc.Card([
                    "Li, J., Hubisz, M., Earlie, E., Duran, M., et al. “Non-cell-autonomous Progression of Cancer Cells From Chromosomal Instability.” ",
                    html.I("Nature"),
                    ", 2023, ",
                    html.A("https://doi.org/10.1038/s41586-023-06464-z", href="https://doi.org/10.1038/s41586-023-06464-z"),
                    "."
                ], body=True, color='secondary')])
    ]


def home_approach_description():
    return [
        dbc.Row([
            dbc.Col(html.Img(className='img-fluid', src=dash.get_asset_url("contacttracing_summary.png"),
                              style={"width": "66%"},
                              alt="ContactTracing Summary"), className='g-0 text-end', align='right'),
            dbc.Col(html.Img(className='img-fluid', src=dash.get_asset_url("contacttracing_summary2.png"),
                              style={"width": "46%"},
                              alt="ContactTracing Summary 2"), className='g-0', align='left'),
        ], className='g-0 justify-content-center'),
        html.P(),
        html.Div([
            "Cells typically send messages to their microenvironment by emitting ligands, which bind to cognate "
            "receptors on the surfaces of target cells; then triggering a change in signaling response in the target cell. "
            "Complimentary ligand-receptor-pairs are well annotated and methods have been developed to infer cell-cell "
            "interactions from single cell expression data based on their mutual expression. ",
            html.I("ContactTracing"),
            " builds upon these existing methods by evaluating ",
            html.I(html.U("the effect")),
            " of such cell-cell interaction on gene expression in the target (receptor-expressing) cell upon ligand "
            "engagement. How does ",
            html.I("ContactTracing"),
            " do this? It makes the simple assumption that not all cells are actively engaged in a putative cell-cell "
            "interaction. By exploiting intrinsic biological variability in (1) receptor expression on target cells and "
            "(2) ligand availability in the tumor microenvironment through a series of nested comparisons, "
            "we infer downstream gene expression changes in the target (receptor-expressing) cell type upon ligand "
            "engagement."
        ])
    ]


def home_dataset_descriptions():
    return [
        html.Div([
            "Foundational work towards benchmarking a validation of this method was performed in isogenic breast cancer "
            "models distinguished by tumor cell-intrinsic rates of chromosome missegregation (a cellular process called "
            "chromosomal instability, or CIN) or its downstream activation of the cGAS-STING antiviral sensing pathway. "
            "Through this, we identified tumor ligands emanating from an endoplasmic reticulum (ER)-stress response as "
            "potential mediators of immune suppression in chromosomally  unstable tumors. Indeed, CIN-induced chronic "
            "STING activation led to rapid interferon-selective desentization and a switch to ER-stress-dependent "
            "transcription. Moreover, inhibition of chronic STING - or key mediators of its unfolded protein response "
            "to ER-stress - suppressed metastasis in syngeneic models of melanoma, breast, and colorectal cancer; "
            "validating this innovative methodology and identtifying a targetable mediator of cancer metastasis! Read "
            "more about our findings ",
            html.A("here", href="https://doi.org/10.1038/s41586-023-06464-z"),
            "."
        ]),
        html.P(),
        html.Div([
            html.H5("Models to study CIN-dependent effects on the tumor microenvironment (TME)"),
            "CIN itself is an ",
            html.I("ongoing process"),
            " of chromosome segregation errors during mitosis. To test whether this ",
            html.I("active process"),
            " itself contributes to disease progression independent of the aneuploidy that arises as a result of it, we "
            "developed experimental tools that enable us to manipulate chromosome missgregation rates in otherwise "
            "isogenic tumor models by altering kinetochore-microtubule attachment stability, (link to papers: ",
            html.A("https://www.nature.com/articles/ncb1809", href="https://www.nature.com/articles/ncb1809"),
            " and ",
            html.A("https://www.nature.com/articles/nature25432", href="https://www.nature.com/articles/nature25432"),
            "). Critically, both CIN",
            html.Sup("high"),
            " and CIN",
            html.Sup("low"),
            " tumors are equally aneuploid, yet differ in their ongoing rates of segregation errors during cell "
            "division. As such, this experimental system enables us to dissect the relative contributions of CIN - as "
            "defined by the ",
            html.I("rate"),
            " of ongoing chromosome missegregation - from aneuploidy, which is a ",
            html.I("state"),
            " of abnormal chromosome numbers. Here, 10X Chromium scRNA-seq data was collected from the 4T1 murine model "
            "of metastatic, triple negative breast cancer, in which we selectively ",
            html.I("dialed up"),
            " (CIN",
            html.Sup("high"),
            ") or ",
            html.I("dialed down"),
            " (CIN",
            html.Sup("low"),
            ") rates of chromosome missegregation in otherwise isogeneic tumor models."
        ]),
        html.P(),
        dbc.Row([
            dbc.Col(html.Img(className='img-fluid', src=dash.get_asset_url("cin_figure.png"),
                              style={"width": "80%"},
                              alt="CIN Diagram"), className='text-center'),
        ]),
        html.P(),
        html.Div([
            html.H5("The impact of CIN-induced STING signaling on the TME"),
            "Previously, we demonstrated ongoing chromosome segregation errors generate rupture-prone micronuclei that "
            "expose genomic double-stranded DNA (dsDNA) to the cytosol, persistently activating the cGAS-STING pathway "
            "to drive metastasis (",
            html.A("https://www.nature.com/articles/nature25432", href="https://www.nature.com/articles/nature25432"),
            "). This begged the question: why doesn't STING activation trigger robust anti-tumor immunity? To quantify "
            "how CIN-induced STING signaling impacts cellular crosstalk in the tumor ecosystem, we knocked down (KD) "
            "STING expression in CIN",
            html.Sup("high"),
            " tumors and applied ",
            html.I("ContactTracing"),
            " (STING",
            html.Sup("WT"),
            " vs. STING",
            html.Sup("KD"),
            " in CIN",
            html.Sup("high"),
            " tumors). All conserved CIN- and STING-dependent interactions were intersected, the resultant set was "
            "defined as the maximum value of each respective aggregate statistic across conditions. This dataset "
            "represents the default interaction set visualized by the Circos and Pairwise Interactions Plots.",
        ]),
        html.P(),
        html.Div(html.Img(className='img-fluid', src=dash.get_asset_url("sting_effects.png"),
                          style={"width": "60%"},
                          alt="STING"), className='text-center'),
        html.P(),
        html.Div([
            html.H5("Data Availability"),
            "The single-cell data used in this study is described in our study (",
            html.A("Li et al, 2023", href="https://doi.org/10.1038/s41586-023-06464-z"),
            ") and can be accessed through the Gene Expression Omnibus (GEO) under accession number ",
            html.A("GSE189856", href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189856"), "."
        ])
    ]


def home_plot_descriptions():
    circos_plot_help_button = dbc.Button(wrap_icon("fa-circle-dot", "More Info"), color="info", n_clicks=0, size='md', id='full-circos-help-button')
    cell_type_plot_help_button = dbc.Button(wrap_icon("fa-arrows-left-right-to-line", "More Info"), color="info", n_clicks=0, size='md', id='full-cell-type-help-button')
    ligand_effect_plot_help_button = dbc.Button(wrap_icon("fa-maximize", "More Info"), color="info", n_clicks=0, size='md', id='full-ligand-effect-help-button')

    circos_plot_offcanvas = dbc.Offcanvas(
        circos_help(),
        id="circos-plot-help-offcanvas",
        title="Circos Plot Information",
        is_open=False,
        placement="end",
        backdrop=True,
        keyboard=True
    )
    cell_type_plot_offcanvas = dbc.Offcanvas(
        interactions_help(),
        id="cell-type-plot-help-offcanvas",
        title="Cell Type Interactions Plot Information",
        is_open=False,
        placement="end",
        backdrop=True,
        keyboard=True
    )
    ligand_effect_plot_offcanvas = dbc.Offcanvas(
        ligand_effects_help(),
        id="ligand-effect-plot-help-offcanvas",
        title="Ligand Effect Plot Information",
        is_open=False,
        placement="end",
        backdrop=True,
        keyboard=True
    )

    return [
        circos_plot_offcanvas, cell_type_plot_offcanvas, ligand_effect_plot_offcanvas,
        html.P(),
        html.Div(["Similarly to many scRNA-seq cellular communication toolkits, ",
                html.I("ContactTracing"),
                " is able to identify specific ligand/receptor pairs that are significantly activated within a "
                "single-cell experiment. However, ",
                html.I("ContactTracing"),
                " can capture more subtle and specific interactions within the TME by comparing experimental "
                "conditions. Additionally, unlike other cellular crosstalk profilers, ",
                html.I("ContactTracing"),
                " can identify how these interactions transcriptionally affect downstream genes. Below is a brief "
                "overview of the figures that allow researchers to interpret the high dimensional output of ",
                html.I("ContactTracing"),
                ":"]),
        html.H5([
            html.Div([
                html.I(html.A("Circos Plot", href="/circos")),
                circos_plot_help_button,
            ], className="d-grid d-md-flex gap-1 justify-content-md-between"),
        ]),
        html.Div(html.Img(className='img-fluid', src=dash.get_asset_url("circos_help.png"), style={"width": "36%"}, alt="Circos Plot"), className='text-center'),
        html.P(),
        html.Div(["The Circos plot summarizes all condition-specific interactions between cells in the TME as identified "
                "by ",
                html.I("ContactTracing"),
                ". The ribbons highlight the strongest ligand/receptor interactions across cell types. The Circos "
                "diagram can be helpful for rapid hypothesis generation from the data."]),
        html.H5([
            html.Div([
                html.I(html.A("Cell Type Interactions Plot", href="/interactions")),
                cell_type_plot_help_button,
            ], className="d-grid d-md-flex gap-1 justify-content-md-between"),
        ]),
        html.Div(html.Img(className='img-fluid', src=dash.get_asset_url("pairwise_help2.png"), style={"width": "34%"}, alt="Cell Type Interactions Plot"), className='text-center'),
        html.P(),
        html.Div(["The pairwise cell type interactions plot highlights the condition-specific signals sent by a donor "
                "cell type (emitting ligands) to a target cell type (expressing receptors). This figure gives a "
                "focused view of cell-cell interactions of interest."]),
        html.H5([
            html.Div([
                html.I(html.A("Downstream Ligand Effects Plot", href="/ligand-effects")),
                ligand_effect_plot_help_button,
            ], className="d-grid d-md-flex gap-1 justify-content-md-between"),
        ]),
        html.Div(html.Img(className='img-fluid', src=dash.get_asset_url("ligand_effects_help3.png"), style={"width": "24%"}, alt="Downstream Ligand Effects Plot"), className='text-center'),
        html.P(),
        html.Div("This figure illustrates the downstream effects of a ligand within the TME. Given a donor cell type and "
               "ligands of interest, it is possible to visualize each of the receptors that are activated in a "
               "CIN-dependent manner and how these receptors induce cascading transcriptional responses throughout "
               "the microenvironment.")
    ]


def home_misc_info():
    contact_info_button = dbc.Button(wrap_icon("fa-envelope", "View Email"), color="info", n_clicks=0, size='lg', id='contact-info-button')
    contact_div = html.Div([], id="contact-info")

    return html.Div([
        html.P(),
        html.Div(["Below are some details regarding the implementation of ",
                html.I("ContactTracing"),
                " and the dashboards presented on this website. For more details, please refer to our paper."]),
        html.P(),
        html.H5([html.I("ContactTracing"), " Information"]),
        html.P(),
        html.Div([html.I("ContactTracing"),
                " uses the ligand/receptor pairing databases from ",
                html.A("CellTalkDB", href="http://tcm.zju.edu.cn/celltalkdb/download.php"),
                " and ",
                html.A("CellphoneDB", href="https://www.cellphonedb.org/database.html"),
                ". In addition, the interaction tests are implemented using the ",
                html.A("MAST framework", href="https://rglab.github.io/MAST/"),
                ", which fits a zero-inflated hurdle model using our model design specifications to the scRNA-seq "
                "expression counts data. "]),
        html.P(),
        html.H5("Website Information"),
        html.P(),
        html.Div(["The website is implemented using the ",
                html.A("Plotly Dash framework", href="https://dash.plotly.com/"),
                ", using the pre-computed ",
                html.I("ContactTracing"),
                " analysis data from our corresponding paper. This website aims to emulate a few of the figures "
                "presented in our work but allows for user interactivity. A slightly modified version of these interactive"
                " figures can be generated by following the jupyter notebook found on the website's ",
                html.A("GitHub page", href="https://github.com/LaughneyLab/ContactTracing-Viz"),
                ". To learn more about the Laughney Lab, check us out at: ",
                html.A("https://laughneylab.com/", href="https://laughneylab.com/"),
                ". If any issues are encountered with this "
                "site, or there are any other concerns, please contact us using the information below:"]),
        html.P(),
        html.H5("Contact Information"),
        html.P(),
        html.Div(contact_info_button),
        html.P(),
        html.Div(contact_div)
    ])


def conditions_def(text):
    return make_tooltip(text, [
        html.I("ContactTracing"),
        " compares the transcriptional profiles of genes between two conditions. In our study, we compared "
        "CINlow vs CINhigh tumors and CINhigh vs STING-KD tumors. Presented are the results for the CINlow vs CINhigh "
        "comparison and the intersection between the CINlow vs CINhigh and CINhigh vs STING-KD comparisons."
    ])


def ligand_log2fc_def(text):
    return make_tooltip(text, [
        "As part of the tests ",
        html.I("ContactTracing"),
        " performs, differential availability of a ligand between conditions is required to identify condition-specific "
        "effects from the TME."
    ])


def diffusion_component_def(text):
    return make_tooltip(text, "The Diffusion component metric depicted is an embedding of each gene's differential expression, "
                              "it is used for ordering genes along the Circos plot.")


def differential_abundance_def(text):
    return make_tooltip(text, "The differential abundance metric depicted represents the relative correlation of each genes expression "
                              "across each cell type's differential abundance between conditions.")


def deg_test_def(text):
    return make_tooltip(text,
                        [
                            html.I("ContactTracing"),
                            " uses MAST to identify whether a downstream gene's expression is significantly differentially"
                            " expressed according to the activation of a given receptor."
                        ])


def interaction_test_def(text):
    return make_tooltip(text,
                        [
                            html.I("ContactTracing"),
                            " uses MAST to test for whether a given gene's expression is significantly dependent on "
                            "the condition-specific activation of a given receptor."
                        ])


def interaction_effects_def(text):
    return make_tooltip(text,
                        "An interaction effect is a gene which has expression that is effected by "
                        "condition-specific activation of a given receptor.")


def circos_help():
    return [
        html.H5(["The ",
                 html.I('ContactTracing'),
                 " Circos plot aggregates all relevant interaction information into a single figure."]),
        html.Hr(),
        html.P(),
        html.Div(html.H6(html.Strong("Filter Options"))),
        html.Ul([
            html.Li([html.P(),
                html.Div([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a receptor"
                " within a given cell of interest; this modulates the threshold for classifying induced expression"
                " change as significant after a Benjamini-Hochberg FDR correction. The default value of 0.25 reflects"
                " what was chosen for evaluating the intersection between CIN and STING-dependent effects. A cutoff"
                " of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Ligand Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific"
                " effects from the Tumor Microenvironment. This cutoff represents the threshold for identifying whether"
                " a ligand's differential expression between conditions significantly differs from 0 after a"
                " Benjamini-Hochberg FDR correction."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Minimum Interaction Effect: As ",
                html.I("ContactTracing"),
                " ContactTracing identifies the condition-specific effects of each receptor's activation by evaluating "
                "all possible genes, this value indicates a filter for the minimum number of condition-specific "
                "activations of genes by each given receptor using the selected Interaction Effect FDR Cutoff. The "
                "interaction effect, in part, determines whether a receptor has interactions depicted in the "
                "\"ribbons\" of the Circos plot."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Minimum Affected Downstream Genes: The number of differentially expressed genes in response to"
                " receptor activation. This value is similar to the Minimum Interaction Effect parameter, but rather"
                " than requiring condition-specific responses, this measures the strength of the Log2FC induced by a"
                " receptor's activation. Note that the FDR-adjusted p-value threshold is the same as the Minimum"
                " Interaction Effect parameter. "
            ])]),
            html.Li([html.P(),
                html.Div([
                "Chord Minimum Ligand abs(Log2FC): Whereas the Ligand Log2FC FDR Cutoff determines whether genes get "
                "included in the plot, this parameter only affects whether a ribbon is drawn to indicate strong "
                "interaction effects within the Circos figure. Increasing this value leads to fewer ribbons drawn "
                "between ligands and receptors, while the opposite is true if decreased."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Genes of Interest: If specified, the chords drawn corresponding to a ligand or receptor listed will "
                "be emphasized (below is an example). Multiple genes can be included by separating names with commas. ",
                html.Em("IMPORTANT: This field is CASE SENSITIVE.")
            ])]),
            html.Li([html.P(),
                html.Div([
                "Interaction Set: We have evaluated both CIN-dependent and CIN & STING-dependent interaction effects. "
                "This toggle lets users instantly see the plot under both conditions. Note that the CIN & STING "
                "effects represent the maximum value between shared interactions across both conditions and require "
                "interaction effects to have the same directionality."
            ])])
        ]),
        html.P(),
        html.Div(html.H6(html.Strong("How to Interpret the Circos Plot"))),
        html.P(),
        html.Div("The final figure is very information-dense, so it can be challenging to interpret initially. "
               "Thus, we will explore the diagram layer by layer. "),
        html.P(),
        html.Div("First, we have sections of the rings grouped by cell type annotation. NOTE: If a cell type name is too "
               "large, it is excluded. However, it is possible to identify the cell type by hovering over each section "
               "of the plot in the browser."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help1.png'), style={"width": '40%', 'align': 'center'},
                          alt="Cell Types", className='img-fluid mx-auto'), className='text-center'),
        html.P(),
        html.Div(["The next layer represents the value of the first Diffusion Component, as calculated by ",
                 html.A("Palantir", href="https://www.nature.com/articles/s41587-019-0068-4"),
               ". This value represents the euclidean space embedding of the differential expression score (the Log2FC "
               "between conditions multiplied by the negative Log10-transformed p-value) of a given gene for a cell "
               "type. The diffusion component allows for a consistent ordering of genes along the Circos rings."]),
        html.Div(html.Img(src=dash.get_asset_url('circos_help2.png'), style={"width": '40%', 'align': 'center'},
                          alt="Diffusion Components", className='img-fluid mx-auto'), className='text-center'),
        html.P(),
        html.Div(["Following Diffusion Component 1, the next ring represents the gene-wise differential abundance of a "
               "gene between conditions. The differential abundance is calculated as the Pearson correlation "
               "coefficient of each gene's imputed expression across cells against the corresponding cell type's "
               "differential abundance across conditions as calculated by ",
                 html.A("MILO", href="https://www.nature.com/articles/s41587-021-01033-z"),
                 ". This value is then normalized to range from -1 to 1."]),
        html.Div(html.Img(src=dash.get_asset_url('circos_help3.png'), style={"width": '40%', 'align': 'center'},
                          alt="Differential Abundance", className='img-fluid mx-auto'), className='text-center'),
        html.P(),
        html.Div("The final ring represents a histogram scaled to the number of significant interaction effects for each "
               "receptor."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help4.png'), style={"width": '40%', 'align': 'center'},
                          alt="Interaction Effects", className='img-fluid mx-auto'), className='text-center'),
        html.P(),
        html.Div("Lastly, all interactions that meet the selected requirements and have the same directionality are "
               "depicted as ribbons connecting each receptor to its associated ligands across cell types. Depending on "
               "the number of ribbons, it may be difficult to read some gene labels. Therefore it is possible to hover "
               "over each ribbon to get precise details. The color of these ribbons reflects a ligand's Log2FC between "
               "conditions, and the thickness of the ribbon is relative to the number of significant interaction "
               "effects for the receptor."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help5.png'), style={"width": '40%', 'align': 'center'},
                          alt="Ribbons", className='img-fluid mx-auto'), className='text-center'),
        html.P(),
        html.Div("Additionally, the ribbons can be filtered to emphasize specific ligands or receptors of interest using "
               "the aforementioned \"Genes of Interest\" field."),
        html.Div(html.Img(src=dash.get_asset_url('circos_help6.png'), style={"width": '40%', 'align': 'center'},
                          alt="Ribbons Highlighted", className='img-fluid mx-auto'), className='text-center')
    ]


def ligand_effects_help():
    return [
        html.H5(["The ligand cascade effects plot visualizes all potential downstream effects from TME-specific"
                 " identified by ",
                 html.I("ContactTracing"),
                 "."]),
        html.Hr(),
        html.P(),
        html.Div(html.H6(html.Strong("Filter Options"))),
        html.Ul([
            html.Li([html.P(),
                html.Div([
                "Emitting Cell Type: The cell type for which specified ligands are emitted."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Emitted Ligands: The ligand(s) to follow in the figure; must be emitted from the selected cell type."
                " Multiple genes can be included by separating names with commas. IMPORTANT: This field is CASE"
                " SENSITIVE."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Network Building Iterations: The number of levels of downstream interactions to include. Note that"
                " typically, the number of downstream genes per level can increase exponentially. So starting with a"
                " low value and increasing as needed is recommended since larger values require additional time to"
                " render."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a "
                " receptor within a given cell of interest; this modulates the threshold for classifying induced"
                " expression change as significant after a Benjamini-Hochberg FDR correction. The default value of 0.25"
                " reflects what was chosen for evaluating the intersection between CIN and STING-dependent effects. A"
                " cutoff of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific"
                " effects from the Tumor Microenvironment. Additionally, as ",
                html.I("ContactTracing"),
                " is able to determine the induced Log2FC of a receptor's activation on a downstream gene, this option"
                " will also be applied to induced Log2FC tests. The cutoff represents the threshold for identifying whether"
                " a gene's differential expression between conditions significantly differs from 0 after a"
                " Benjamini-Hochberg FDR correction."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Minimum Expression: Allows for Ligands/Receptors to be filtered according to a minimum expression"
                " value. This value, however, does not represent counts. Instead, \"expression\" refers to the fraction"
                " of the cells from a given cell type with non-zero counts of a particular gene. Therefore, expression"
                " values range from 0-1, where zero means that no cells of a cell type express the gene, and one means"
                " that all cells of a cell type express the gene."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Minimum abs(Log2FC): While the Log2FC FDR Cutoff option filters interactions according"
                " to whether a gene's differential expression between conditions significantly differs from zero,"
                " this filter additionally allows for further refinement by requiring a minimum absolute Log2FC. The"
                " default value reflects the cutoff used for the Circos plot."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Interaction Set: While we have evaluated both CIN-dependent and CIN & STING-dependent interaction"
                " effects, the CIN & STING intersection data set does not have statistical power for examining"
                " individually induced downstream effects of interactions. Therefore, only CIN-dependent interaction"
                " effects are enabled."
            ])])
        ]),
        html.P(),
        html.Div(html.H6(html.Strong("How to Interpret the Plot"))),
        html.P(),
        html.Div("This figure illustrates the potential cascading effects of receptor activation in the tumor"
               " microenvironment. It does so by generating a directed network of gene interactions by following the"
               " steps listed below:"),
        html.Ol([
            html.Li([html.P(),
                html.Div([
                    "Starting from a given cell type and a set of ligands, identify all the corresponding receptors that"
                    " pass the provided filters and have a significant condition-specific interaction effect on downstream"
                    " signaling, then draw an arrow connecting the ligands to the receptors. Finally, add 1 to the number"
                    " of network building iterations."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help1.png"), style={"width": '30%', 'align': 'center'}, alt='Step 1', className='img-fluid mx-auto'), className='text-center')
            ]),
            html.Li([html.P(),
                html.Div([
                    "Identify significantly differentially regulated ligands for each added receptor in response to the"
                    " receptor's activation. Additionally, identify further receptors with a significant interaction effect"
                    " and a ligand available in the microenvironment. Finally, add 1 to the number of network building"
                    " iterations."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help2.png"), style={"width": '25%', 'align': 'center'}, alt='Step 2', className='img-fluid mx-auto'), className='text-center')
            ]),
            html.Li([html.P(),
                html.Div([
                    "Repeat steps 1 and 2 until the number of network building iterations equals what was specified by"
                    " the user."
                ]),
                html.Div(html.Img(src=dash.get_asset_url("ligand_effects_help3.png"), style={"width": '25%', 'align': 'center'}, alt='Step 3', className='img-fluid mx-auto'), className='text-center')
            ])
        ]),
        html.P(),
        html.Div(["As depicted above, each of these iterations is represented by a column of ligands and receptors. In"
                " addition, these nodes have arrows connecting to their corresponding receptors and ligands. If the"
                " arrow is connecting from a receptor to a ligand, the color and arrow thickness represents the Log2FC"
                " induced by receptor activation. However, if the arrow is connecting from a ligand to a receptor, the"
                " color is equivalent to a Log2FC of 0 as ligands are not expected to regulate receptor expression"
                " differentially."]),
        html.P(),
        html.Div(["Additionally, the color of nodes represents the cell type in which each gene is expressed; the shape"
                " of the node represents whether the gene is a ligand, receptor, or both; and the size of each node is"
                " relative to the total induced Log2FC from receptor activation to ligands."])
    ]


def interactions_help():
    return [
        html.H5(["The pairwise interactions plot highlights TME-specific interactions identified by ",
                 html.I("ContactTracing"),
                 " between cells of interest."]),
        html.Hr(),
        html.P(),
        html.Div(html.H6(html.Strong('Filter Options'))),
        html.Ul([
            html.Li([html.P(),
                html.Div([
                "Interaction Effect FDR Cutoff: ",
                html.I("ContactTracing"),
                " identifies the downstream genes that have induced expression shifts from the activation of a "
                "receptor within a given cell of interest; this modulates the threshold for classifying induced "
                "expression change as significant after a Benjamini-Hochberg FDR correction. The default value of "
                "0.25 reflects what was chosen for evaluating the intersection between CIN and STING-dependent "
                "effects. A cutoff of 0.05 might be more sensible when evaluating CIN-dependent effects on its own."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Ligand Log2FC FDR Cutoff: ",
                html.I("ContactTracing"),
                " requires differential availability of a ligand across conditions to identify condition-specific "
                "effects from the Tumor Microenvironment. This cutoff represents the threshold for identifying "
                "whether a ligand's differential expression between conditions significantly differs from 0 after a "
                "Benjamini-Hochberg FDR correction."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Minimum Interaction Effect: As ",
                html.I("ContactTracing"),
                " identifies the condition-specific effects of each receptor's activation by evaluating all possible "
                "genes, this value indicates a filter for the minimum number of condition-specific activations of "
                "genes by each given receptor using the selected Interaction Effect FDR Cutoff."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Minimum Expression: Allows for Ligands/Receptors to be filtered according to a minimum expression "
                "value. This value, however, does not represent counts. Instead, \"expression\" refers to the "
                "fraction of the cells from a given cell type with non-zero counts of a particular gene. Therefore, "
                "expression values range from 0-1, where zero means that no cells of a cell type express the gene, "
                "and one means that all cells of a cell type express the gene."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Minimum Ligand abs(Log2FC): While the Ligand Log2FC FDR Cutoff option filters interactions according "
                "to whether a ligand's differential expression between conditions significantly differs from zero, "
                "this filter additionally allows for further refinement by requiring a minimum absolute Log2FC. The "
                "default value reflects the cutoff used for the Circos plot."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Interaction Directionality: By default, the figure generated only depicts pairwise interactions "
                "between cell types from left to right. In other words, for every cell type depicted, interactions "
                "shown depict ligands emitted from the cell type directly to the left of a given cell type. "
                "Bidirectional interactions can also be enabled, which can visualize the two-way cross-talk of cells."
            ])]),
            html.Li([html.P(),
                html.Div([
                "First Cell Type: The first cell type of interest to depict interactions between."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Second Cell Type: The second cell type of interest to depict interactions between."
            ])]),
            html.Li([html.P(),
                html.Div([
                "Third Cell Type: The third cell type to include in the figure; this is optional. If not specified, "
                "only two cell types will be shown.",
            ])]),
            html.Li([html.P(),
                html.Div([
                "Interaction Set: We have evaluated both CIN-dependent and CIN & STING-dependent interaction effects. "
                "This toggle lets users instantly see the plot under both conditions. Note that the CIN & STING "
                "effects represent the maximum value between shared interactions across both conditions and require "
                "interaction effects to have the same directionality."
            ])])
        ]),
        html.P(),
        html.Div(html.H6(html.Strong('How to Interpret the Plot'))),
        html.P(),
        html.Div("Selected cell types will be represented as columns in the figure according to the order selected "
               "(ex. first cell type is the leftmost column), with selected genes listed within it. These genes may "
               "be either a ligand (circle), receptor (square), or a gene that is both a ligand and receptor (diamond "
               "with dot). Additionally, these genes are colored according to the proportion of cells of the given "
               "cell type expressing each gene. Note: when many genes are selected, it is possible for labels to "
               "overlap, in which case the user can toggle the label of each gene by clicking on the representative "
               "node. Columns of genes are illustrated below:"),
        html.Div(html.Img(src=dash.get_asset_url('pairwise_help1.png'), style={"width": '20%', 'align': 'center'}, alt='Cell Type Columns', className='img-fluid mx-auto'), className='text-center'),
        html.P(),
        html.Div("To illustrate interactions between ligands and receptors, arrows are drawn between pairs that have "
               "interactions that meet the user-defined filters. The colors of these arrows indicate the Log2FC "
               "strength and directionality of the ligand between conditions (red corresponds to up-regulation, and "
               "blue corresponds to down-regulation). Additionally, each arrow has a thickness relative to the number "
               "of significant interaction effects in the receptor. Below is an example of the arrows in a "
               "unidirectional interactions plot. However, users can allow interactions to start at either cell type "
               "by selecting the \"bidirectional\" Interaction Directionality setting."),
        html.Div(html.Img(src=dash.get_asset_url('pairwise_help2.png'), style={"width": '20%', 'align': 'center'}, alt='Arrows', className='img-fluid mx-auto'), className='text-center')
    ]
