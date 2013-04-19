#!/usr/bin/ruby

require File.dirname(__FILE__) + "/../common/libgraphload.rb"


NNODE = ARGV[0].to_i
NPATH = ARGV[1].to_i
SEEDS = ARGV[2].to_i
OUTPUTFILE = ARGV[3].to_s

wout = open(OUTPUTFILE, "w")
abort "open file ${OUTPUTFILE} error" unless wout


#read topology
C = topology(NNODE, SEEDS)
D = traffic(NNODE, NPATH)


#get we with edge list style
$e = []
NNODE.times { |i| NNODE.times { |j| next if i == j
    next if C[i][j] == 0
    $e << [i,j]
} }
NEDGE = $e.length

wout.puts "#define N #{NNODE}"
wout.puts "#define E #{NEDGE}"
wout.puts "#define H #{NPATH}"
wout.puts "#define S #{SEEDS}"



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






















