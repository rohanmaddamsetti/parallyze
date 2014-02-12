#!/usr/bin/python

#usage: python parallyze.py config.txt [genomediff1.gd genomediff2.gd ... n : referencegenome.gb(k) {genbank}]

#find config.txt
#find .gd files, save in list
#try to break

import argparse

def main():
	parser = argparse.ArgumentParser()
#	parser.parse_args()
	parser.add_argument(dest='config', help="First entry should be a config.txt file. File should be in working folder.")
	parser.add_argument('-r', dest='reference', default="NONE")
	gdfiles=[]
      #  for i in cmd/argument/args:
       #         if i.endswith(".gd"):
        #                gdfiles.append(i)
	args = parser.parse_args()
	
	print "\nHello World!\n"
	if args.config != 'config.txt':
		print "Error: please ensure file is titled config.txt"
#	elif 'config'!=a text file in right format
#		print "Error: please ensure file is a .txt file"
#	elif 'config' != contents of file properly formatted/layed out
#		print "Error: please ensure contents of file are appropriate"
	if args.reference:
		print "The reference genome is", args.reference
main()
