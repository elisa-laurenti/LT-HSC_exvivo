suppressMessages(library(gdata))
suppressMessages(library(googleVis))
suppressMessages(library(ggsankey))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(htmlwidgets))
suppressMessages(library(htmltools))

sankey_val = read.csv('dataframe_for_sankey_plotting_BASiCS_DV_N_geneTrends.csv',row.names=1)

colors_link <- c('blue', 'cornflowerblue', 'green', 'red')
colors_link_array <- paste0("[", paste0("'", colors_link,"'", collapse = ','), "]")
colors_node <- c('black', 'black', 'black', 'black')
colors_node_array <- paste0("[", paste0("'", colors_node,"'", collapse = ','), "]")
opts <- paste0("{
        link: { colorMode: 'source',
                colors: ", colors_link_array ," },
        node: { colors: ", colors_node_array ," },
        iterations: 0
      }" 
               )

sk2 <- gvisSankey(sankey_val, from="Source", to="Target", weight="Value",
                 options=list(
                   sankey=opts))
