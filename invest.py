#!/usr/bin/env python
#program to compute value of investment 10 years in future

def main():
	print
	print("This program calculates future value of a 10-year investment")
	print
	principal=input("What is your principal? : ")
	apr=input("What is the annual percentage rate? : ")
	for i in range(10):
		principal=principal*(1+apr)
	print
	print "The value of your principal after 10 years is :", principal
	print
main()
