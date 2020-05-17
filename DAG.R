library(DiagrammeR)
library(rsvg)
library(DiagrammeRsvg)

#embed windows fonts
#Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.52/bin/gswin64c.exe")


#Black point DAG showing the interaction terms
graph <- grViz("

digraph boxes_and_circles {

# add node statements
node [shape  = circle,
      fillcolor = grey,
      penwidth = 1.5,
      style = filled]

NEC; NEP; pH; Day_Night[label = 'Day/Night']; Tide; NN; Season; Temperature; Day_Night_Tide[label = 'Day/Night x \nTide'] ;
Day_Night_NN[label = 'Day/Night x \nNN']; Tide_NN[label = 'Tide x NN']; Season_Temperature[label = 'Season x \nTemperature']; 
Day_Night_Tide_NN[label = 'Day/Night x \nTide x \nNN']; SGD[label = '%SGD']; 
Day_Night_SGD[label = 'Day/Night x \n%SGD']; Day_Night_Season[label = 'Day/Night x \nSeason'];
SGD_Season[label = '%SGD x \nSeason']; Day_Night_SGD_Season[label = 'Day/Night x \n%SGD x \nSeason'];

# add edge statements
SGD -> Temperature[color = red, penwidth = 0.24]
SGD -> NN [color = blue, penwidth = 1.1]
SGD -> pH [color = blue, penwidth = 0.17]
Season -> Temperature[color = red, penwidth = 1.4]
Day_Night -> Temperature[color = red, penwidth = 0.76]
Day_Night_SGD -> Temperature[color = blue, penwidth = 0.3]
Day_Night_Season -> Temperature[color = red, penwidth = 0.26]
SGD_Season-> Temperature[color = grey, penwidth = 0.1]
Day_Night_SGD_Season-> Temperature[color = grey, penwidth = 0.15]
NN -> NEP [color = blue, penwidth = 0.41]
Season -> NEP[color = blue, penwidth = 2]
Tide -> NEP[color = red, penwidth = 0.57]
Day_Night -> NEP[color = red, penwidth = 0.74]
Temperature -> NEP[color = blue, penwidth = 1]
Day_Night_Tide -> NEP[color = blue, penwidth = 0.57]
Day_Night_NN -> NEP[color = red, penwidth = 0.63]
Tide_NN -> NEP[color = red, penwidth = 0.74]
Season_Temperature -> NEP[color = grey, penwidth = 0.19]
Day_Night_Tide_NN -> NEP[color = blue, penwidth = 0.7]
NEP-> pH[color = blue, penwidth = 0.93]
Temperature -> NEC[color = grey, penwidth = 0.1]
pH -> NEC[color = blue, penwidth = 0.59]

# add a graph statement
graph [nodesep = 0.1]

}
")

# DAG showing marginal effects
graph_BlackPoint <- grViz("

digraph boxes_and_circles {

# add node statements
node [shape  = circle,
      fillcolor = grey,
      penwidth = 1.5,
      style = filled]

NEC; NEP; pH; Day_Night[label = 'Day/Night']; Tide; NN; Season; Temperature; SGD[label = '%SGD'];

# add edge statements
SGD -> Temperature[color = red, penwidth = 0.24, label = 'Fall Day']
SGD -> Temperature[color = red, penwidth = 0.16, label = 'Spring Day']
SGD -> Temperature[color = red, penwidth = 0.06, label = 'Fall Night']
SGD -> Temperature[color = red, penwidth = 0.01, label = 'Spring Night']
Season -> Temperature[color = red, penwidth = 1.4]
Day_Night -> Temperature[color = red, penwidth = 0.76]

SGD -> NN [color = blue, penwidth = 1.1]

SGD -> pH [color = blue, penwidth = 0.17]

Season -> NEP[color = blue, penwidth = 2]
Tide -> NEP[color = red, penwidth = 0.57]
Day_Night -> NEP[color = red, penwidth = 0.74]

NN -> NEP [color = blue, penwidth = 0.42, label = 'Day High']
NN -> NEP [color = red, penwidth = 0.33, label = 'Day Low']
NN -> NEP [color = red, penwidth = 0.23, label = 'Night High']
NN -> NEP [color = red, penwidth = 0.26, label = 'Night Low']

Temperature -> NEP[color = blue, penwidth = 1, label = 'Fall']
Temperature -> NEP[color = blue, penwidth = 1.2, label = 'Spring']

NEP-> pH[color = blue, penwidth = 0.93]

Temperature -> NEC[color = grey, penwidth = 0.1]
pH -> NEC[color = blue, penwidth = 0.59]

# add a graph statement
graph [nodesep = 0.1]

}
")

#export
graph_BlackPoint %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/BlackPointDAG.pdf")

#Wailupe
graph_Wailupe <- grViz("

digraph boxes_and_circles {

# add node statements
node [shape  = circle,
      fillcolor = grey,
      penwidth = 1.5,
      style = filled]

NEC; NEP; pH; Day_Night[label = 'Day/Night']; Tide; NN; Season; Temperature; SGD[label = '%SGD'];

# add edge statements
SGD -> Temperature[color = red, penwidth = 0.51, label = 'Fall Day']
SGD -> Temperature[color = red, penwidth = 0.41, label = 'Spring Day']
SGD -> Temperature[color = red, penwidth = 0.35, label = 'Fall Night']
SGD -> Temperature[color = red, penwidth = 0.23, label = 'Spring Night']
Season -> Temperature[color = red, penwidth = 1.3]
Day_Night -> Temperature[color = red, penwidth = 0.61]

SGD -> NN [color = blue, penwidth = 0.67]

SGD -> pH [color = red, penwidth = 0.11]

Season -> NEP[color = red, penwidth = 0.58]
Tide -> NEP[color = grey, penwidth = 0.35]
Day_Night -> NEP[color = red, penwidth = 1.4]

NN -> NEP [color = red, penwidth = 0.07, label = 'Day High']
NN -> NEP [color = blue, penwidth = 0.33, label = 'Day Low']
NN -> NEP [color = red, penwidth = 0.05, label = 'Night High']
NN -> NEP [color = blue, penwidth = 0.17, label = 'Night Low']

Temperature -> NEP[color = blue, penwidth = 0.26, label = 'Fall']
Temperature -> NEP[color = blue, penwidth = 0.18, label = 'Spring']

NEP-> pH[color = blue, penwidth = 0.95]

Temperature -> NEC[color = blue, penwidth = 0.22]
pH -> NEC[color = blue, penwidth = 0.24]

# add a graph statement
graph [nodesep = 0.1]

}
")

graph_Wailupe %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/WailupeDAG.pdf")


