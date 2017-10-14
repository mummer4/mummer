#! /usr/bin/env ruby

require "mummer"

ref = File.read(ARGV[0])
qry = File.read(ARGV[1])

o = Mummer::Options.new.minmatch(10).mincluster(10)
aligns = Mummer::align_sequences(ref, qry, o)

aligns.each { |a|
  puts("#{a.sA} #{a.eA} #{a.sB} #{a.eB} #{a.Errors} #{a.SimErrors} #{a.NonAlphas}")
  puts(a.delta.join("\n")) unless a.delta.empty?
  puts("0")
}
