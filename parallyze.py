#!/usr/bin/python

#usage: python parallyze.py config.txt [genomediff1.gd genomediff2.gd ... n : referencegenome.gb(k) {genbank}]

import argparse
#find config.txt
#find .gd files, save in list
#try to break

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', dest='reference', default="NONE")
	parser.add_argument(dest='config')
	gdfiles=[]
        for i in cmd:
                if i.endswith(".gd")
                        gdfiles.append(i)
	args = parser.parse_args()
	
	print "Hello World!"
	if args.reference:
		print "The reference genome is", args.reference
main()
