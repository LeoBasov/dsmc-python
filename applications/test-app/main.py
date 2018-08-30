import sys
sys.path.append('../../src')

import dsmc
from plasma import generator
from IO.input import read_xml

def print_header():
	print(40*"=")
	print('DSMC python test')
	print(40*"-")

def print_footer():
	print(40*"=")

def main():
	print_header()
	
	read_xml(sys.argv[1])

	print_footer()

main()