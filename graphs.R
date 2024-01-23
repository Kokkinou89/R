install.packages("igraph") 
library(igraph)
setwd("C:/")

#Graph 1

rm(list = ls())

ed_list = read.table('C:/Users/user/Desktop/–·Ò‹‰ÔÛÁ/ œ  …Õœ’ –¡Õ¡√…Ÿ‘¡/≈Ò„·ÛÈ·1/molecular_data.txt')
ed_list = as.matrix(ed_list)

vertex_name = read.table('C:/Users/user/Desktop/–·Ò‹‰ÔÛÁ/ œ  …Õœ’ –¡Õ¡√…Ÿ‘¡/≈Ò„·ÛÈ·1/molecular_names.txt')
vertex_name = as.matrix(vertex_name)

net = graph.edgelist(ed_list,directed=TRUE)

plot(net,vertex.label.font = 2,
     vertex.size = 5,
     vertex.color = rainbow(10, 0.8, 0.8, alpha = 0.7),
     vertex.label = vertex_name,
     vertex.label.color = "black",
     vertex.frame.color = 'grey',
     vertex.label.cex = 0.6,
     vertex.label.degree = -pi/2,
     edge.arrow.size = 0.25,
     main='Steroid Hormone Biosynthesis for Mus Musculus (mmu00140)',
     asp = 0, margin = -0.1)

##Degrees

degree_net = as.data.frame(degree(net))

degree_data_net = cbind(vertex_name,
                        degree_net)

colnames(degree_data_net) = c('Molecular Name',
                              'Degree')

degree_data_net

##Shortest path

new_vertex_name = as.data.frame(cbind(V(net),vertex_name))
colnames(new_vertex_name) = c('ID','Molecular Name')

Shortest_path_net = shortest.paths(net, as.character(new_vertex_name[2,1]),
                                   as.character(new_vertex_name[7,1]))

colnames(Shortest_path_net) = new_vertex_name[2,2]

rownames(Shortest_path_net) = new_vertex_name[7,2]

Shortest_path_net

#Metrics

##Betweenness

betweenness_net = as.data.frame(betweenness(net))

betweenness_net = cbind(vertex_name,
                        betweenness_net)

colnames(betweenness_net) = c('Molecular Name',
                              'Betweenness')

betweenness_net[order(betweenness_net$Betweenness,decreasing = T),]

##Closeness

closeness_net = as.data.frame(closeness(net))

closeness_net = cbind(vertex_name,
                      closeness_net)

colnames(closeness_net) = c('Molecular Name',
                            'Clossenness')

closeness_net[order(closeness_net$Clossenness,decreasing = T),]

##Alpha

alpha_net = as.data.frame(alpha.centrality(net))

alpha_net = cbind(vertex_name,
                  alpha_net)

colnames(alpha_net) = c('Molecular Name',
                        'Alpha')

alpha_net[order(alpha_net$Alpha,decreasing = T),]

##K authority

K_authority_net = as.data.frame(authority.score(net))

vertex_name = as.data.frame(vertex_name)

K_authority_net = cbind(vertex_name,
                        data.matrix(K_authority_net[,1]))

colnames(K_authority_net) = c('Molecular Name', "Kleinberg's authority")

K_authority_net[order(K_authority_net$`Kleinberg's authority`,decreasing = T),]


##K Hub

K_hub_net = as.data.frame(hub.score(net))

K_hub_net = cbind(vertex_name,
                  data.matrix(K_hub_net[,1])) 

colnames(K_hub_net) = c('Molecular Name', "Kleinberg's hub")

K_hub_net[order(K_hub_net$`Kleinberg's hub`,decreasing = T),]

#Graph 2

rm(list = ls())

e_l = read.table('C:/Users/user/Desktop/–·Ò‹‰ÔÛÁ/ œ  …Õœ’ –¡Õ¡√…Ÿ‘¡/≈Ò„·ÛÈ·1/brain_data.txt')
adj_matrix = as.matrix(e_l)
colnames(adj_matrix) = c(1:ncol(adj_matrix))

v_n = read.table('C:/Users/user/Desktop/–·Ò‹‰ÔÛÁ/ œ  …Õœ’ –¡Õ¡√…Ÿ‘¡/≈Ò„·ÛÈ·1/brain_label.txt')
v_n = as.matrix(v_n)

g = graph.adjacency(adj_matrix,mode = c("directed"))

plot(g, vertex.label.font = 2,
     vertex.size = 6,
     vertex.color = rainbow(10, 0.8, 0.8, alpha=0.8),
     vertex.label = v_n,
     vertex.label.color = "black",
     vertex.frame.color='grey',
     vertex.label.cex = 0.7,
     vertex.label.degree = -pi/2,
     edge.arrow.size = 0.25,
     main='Brain Topology Network During Mental Rotation Task',
     asp = 0, margin = 0)


##Degrees

degree_g = as.data.frame(degree(g))

degree_data_g = cbind(v_n,
                      degree_g)

colnames(degree_data_g) = c('Brain Label','Degree')

degree_data_g

##Shortest path

new_v_n = as.data.frame(cbind(V(g)$name,v_n))

Shortest_path_g = shortest.paths(g, as.character(new_v_n[45,1]),
                                 as.character(new_v_n[32,1]))

colnames(Shortest_path_g) = new_v_n[45,2]

rownames(Shortest_path_g) = new_v_n[32,2]

Shortest_path_g

#Metrics

##Betweenness

betweenness_g = as.data.frame(betweenness(g))

betweenness_g = cbind(v_n, betweenness_g)

colnames(betweenness_g) = c('Brain Label', 'Betweenness')

betweenness_g[order(betweenness_g$Betweenness,decreasing = T),]

##Closeness

closeness_g = as.data.frame(closeness(g))

closeness_g = cbind(v_n,
                    closeness_g)

colnames(closeness_g) = c('Brain Label', 'Clossenness')

closeness_g[order(closeness_g$Clossenness,decreasing = T),]

##Alpha

alpha_g = as.data.frame(alpha.centrality(g))

alpha_g = cbind(v_n,alpha_g)

colnames(alpha_g) = c('Brain Label','Alpha')

alpha_g[order(alpha_g$Alpha,decreasing = T),]

##K authority

K_authority_g = as.data.frame(authority.score(g))

K_authority_g = cbind.data.frame(v_n, K_authority_g[,1])

colnames(K_authority_g) = c('Brain Label', "Kleinberg's authority")

K_authority_g[order(K_authority_g$`Kleinberg's authority`,decreasing = T),]


##K Hub

K_hub_g = as.data.frame(hub.score(g))

K_hub_g = cbind.data.frame(v_n,K_hub_g[,1]) 

colnames(K_hub_g) = c('Brain Label', "Kleinberg's hub")

K_hub_g[order(K_hub_g$`Kleinberg's hub`,decreasing = T),]