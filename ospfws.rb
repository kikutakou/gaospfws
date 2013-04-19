#!/usr/bin/ruby

require File.dirname(__FILE__) + "/../common/libecmp.rb"
require File.dirname(__FILE__) + "/../common/resultshandling/cplex_readresult.rb"
require File.dirname(__FILE__) + "/../common/libubigraph.rb"
require File.dirname(__FILE__) + "/../common/libroute.rb"



FILE = ARGV[0].to_s
if File.exist?(FILE)
	then
	NNODE,NPATH,SEEDS = FILE.split(/[^0-9]+/).reject!{|x| x.to_i < 1 }.map!{ |x| x.to_i }
	else
	NNODE,NPATH,SEEDS = ARGV[0].to_i, ARGV[1].to_i, ARGV[2].to_i
end

abort "usage :  #{__FILE__}  [outputFile] / NNODE NPATH SEEDS" if NNODE < 3 or NPATH < 1 or !defined?(SEEDS)
DEBUG = ARGV[3].to_i unless ARGV[3].nil?


#read topology
C = topology(NNODE, SEEDS)
D = traffic(NNODE, NPATH)
w = Array.new(NNODE, nil) ; w.map! { |x| Array.new(NNODE, nil) }
$e = [] ; NNODE.times { |i| NNODE.times { |j| next unless C[i][j] > 0 ; $e << [i,j] } }

we = [   6,   3,   4,   5,   3,   3,   6,  18,   4,   3,   2,  16,  15,  11,  19,   8,   5,   9,  10,  13,   9,  16,  10,  14,   3,   7,   8,  17,  10,  17,   4,  14,  15,  20,  12,  18,  16,  14]
abort "we definition worng" unless we.length == $e.length

if File.exist?(FILE)
	then
	##read results
	_,_,_,ww,_,_ = cplexread(FILE)
	ww.map { |x, v| w[x]=Array.new(NNODE, nil) ; v.map { |y, v| w[x][y] = v } }
	$e.each_with_index{ |x,i| we[i] = w[x[0]][x[1]] }
	else
	$e.each_with_index{ |x,i| w[x[0]][x[1]] = we[i] }
end



#puts "weight" if defined?(DEBUG)
#w.map{ |x| p x } if defined?(DEBUG)
#puts "edge"
#printedge($e, we)

#---- DIJKSTRA
#n,d = ecmpdijkstra(w)
#puts "ecmp dijkstra : next table n"
#n.map{ |x| p x }

#---- WARSHALL_FLOYD (use -1 to show local edge)
p,d = ecmpwarshallfloyd(w)
#puts "ecmp warshall-floyd : pred table p" if defined?(DEBUG)
#p.map{ |x| p x } if defined?(DEBUG)

#n = Array.new( NNODE, nil )
#NNODE.times{ |i|
#    n[i] = Array.new( NNODE, nil )
#    NNODE.times{ |j|
#        n[i][j] = neighborNode_byPredNode(i, j, p)
#        n[i][j].uniq!
#        n[i][j].sort!
#    }
#}



done = Array.new( NNODE, nil )
NNODE.times{ |i| done[i] = Array.new( NNODE, nil ) }

NNODE.times{ |i|
    NNODE.times{ |j|
        follow_pred(i, j, p, w, done)
        
    }
}


#puts "ecmp warshall-floyd : next table n"
#n.map{ |x| p x }

#puts "ecmp warshall-floyd : next table p" if defined?(DEBUG)
#p.map{ |x| p x } if defined?(DEBUG)
#done.map{ |x| p x }






#---- WARSHALL_FLOYD_EDGE
#pe,d = ecmpwarshallfloyd_edge(we, $e)
#puts "ecmp warshall-floyd edge : get predEdgeTable pe ($e=#{$e.length})"
#pe.map{ |x| p x }

#ne = Array.new( pe.length, nil )
#ne.length.times{ |i| ne[i] = Array.new( pe.length, nil )
#    ne[i].length.times{ |j|
#        ne[i][j] = neighborEdge_byPredEdge(i, j, pe, $e)
#        ne[i][j].uniq!
#    }
#}


#ne = getNextEdgeTable_fromPredEdgeTable(pe, $e)
#puts "ecmp warshall-floyd edge : get nextEdgeTable ne ($e=#{$e.length})"
#ne.map{ |x| p x }


#--------- print

#puts "edge"
#printedge($e, we)

#puts "route"
#printroute(n)



#-------  calc L from Pred, Weight, Demand and Capacity

NEDGE = $e.length
$fe = Array.new(NEDGE, 0.0)
$f = Array.new( NNODE, nil )
NNODE.times{ |i| $f[i] = Array.new( NNODE, nil ) }

$e.each{ |i,j| $f[i][j] = 0.0 }
#$f.map{ |x| p x }

#--- alloc flow by tracing n
$e_rev = Hash2.new { |hash,key| hash[key] = Hash.new {} }
$e.each_with_index{ |x,k| $e_rev[x[0]][x[1]] = k }
#p e_rev

def assignFlow_byNextNode(h, d, n, flow, t=1)
    
    t *= n[h][d].length       #no of branches
    
    n[h][d].each{ |j|
        #puts "#{h} - #{d} via #{j} flow #{flow} rate #{t}"
        $f[h][j] += flow/t
        assignFlow_byNextNode(j, d, n, flow, t) unless j == d
    }
end

def assignFlow_byNextEdge(h, d, ne, flow, t=1)
    
    t *= ne[h][d].length       #no of branches
    
    ne[h][d].each{ |k|
        j = $e[k][1]
        #puts "#{$e[k][0]} - #{$e[k][1]} += #{flow} / #{t}"
        $fe[k] += flow/t
        assignFlow_byNextEdge(j, d, ne, flow, t) unless j == d
    }
end

NNODE.times{ |i| NNODE.times{ |j| next if i==j
    next if D[i][j] == 0    #don't execute if no demand
    assignFlow_byNextNode(i, j, p, D[i][j].to_f)
    #assignFlow_byNextEdge(i, j, ne, D[i][j].to_f)
} }

puts "assign flow : f" if defined?(DEBUG)
$f.map{ |x| p x } if defined?(DEBUG)



l = []

$e.each_with_index{ |x,k|
    #puts "edge(#{k}) #{x[0]} - #{x[1]}"
    #puts "#{$f[x[0]][x[1]]} / #{C[x[0]][x[1]]}"
    $f[x[0]][x[1]] /= C[x[0]][x[1]]
    l << $f[x[0]][x[1]]
    #puts "#{$fe[k]} / #{C[x[0]][x[1]]}"
    #l << $fe[k] / C[x[0]][x[1]]
}

puts "assign flow : f" if defined?(DEBUG)
$f.map{ |x| p x } if defined?(DEBUG)

printf("L = %5f\n", l.max)


#puts "p table"
#p.map{ |x| p x }
#puts "d table"
#d.map{ |x| p x }





exit 0




#output to C headder file

HEADDERFILE = File.dirname(__FILE__) + "/w#{NNODE}_#{NPATH}_#{SEEDS}.h"
wout = open(HEADDERFILE, "w")

wout.puts "#define N #{NNODE}"
wout.puts "#define E #{NEDGE}"


#w[i][j]
tmp = []
NNODE.times{ |i| w[i][i] = "0" }
w.map!{ |x| x = x.map!{ |y| y || y="0" } }        #nil expression
w.each{ |x| tmp << "{ " + x.join(", ") + " }" }
wout.puts "const int weight[#{NNODE}][#{NNODE}] = { \n " + tmp.join(", \n ") + "\n};"

#we[k]
wout.puts "const int weightEdge[#{NEDGE}] = { " + we.join(", ") + " };\n"

#C[i][j]
tmp = []
C.each{ |x| tmp << "{ " + x.join(", ") + " }" }
wout.puts "const int capacity[#{NNODE}][#{NNODE}] = { \n " + tmp.join(", \n ") + "\n};"

#D[i][j]
tmp = []
D.each{ |x| tmp << "{ " + x.join(", ") + " }" }
wout.puts "const int demand[#{NNODE}][#{NNODE}] = { \n " + tmp.join(", \n ") + "\n};"

#e[k]
tmp = []
$e.each{ |x| tmp << "{ " + x.join(", ") + " }" }
wout.puts "const int edge[#{NEDGE}][2] = { " + tmp.join(", ") + " };\n"


#puts HEADDERFILE




















