#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim:filetype=python:
#
# Written by Jonas Juselius <jonas.juselius@uit.no> 
# University of Tromsø, 2008
#

import sys
sys.path.append("../Python")
import os, re, optparse
from copy import deepcopy
import getkw

def main():
	print "Starting " + sys.argv[0]
	global topsect
	usage="usage: %prog [options] [config] inpfile"
	cmdln=optparse.OptionParser(usage=usage)
	cmdln.add_option('-V','--version', action='store', dest='version',
			help='print version')
	cmdln.add_option('-v','--verbose', action='store', dest='verbose',
			help='be verbose')
	cmdln.add_option('-d','--debug', action='store', dest='debug',
			help='debug level')
	(opts, args)=cmdln.parse_args()

	if (len(args) == 0):
		inpfil=None
	elif (len(args) == 1):
		inpfil=args[0]
	else:
		cmdln.error('incorrect number of files')
		sys.exit(0)
		
	top=getkw.Section('top', callback=verify_top)
	top.set_status(True)
	top.add_kw('int',			'INT', 0)
	top.add_kw('bool',			'BOOL', True)
	top.add_kw('double',		'DBL', 0.0)
	top.add_kw('int_array',   	'INT_ARRAY')    # array of any lenght
	top.add_kw('double_vector',	'DBL_ARRAY', 3) # array of 3 real
	top.add_kw('string',		'STR')
	top.add_kw('data',			'DATA')

	sect1=getkw.Section('sect1')
	sect1.add_kw('foo',	'INT', 42)
	sect1.add_kw('bar',	'DBL', 74.0e-6)
	sect1.add_kw('oof',	'BOOL', False)
	top.add_sect(sect1)

	subsect=getkw.Section('hello')
	subsect.add_kw('hello',	'STR', "Goodbye")
	subsect.add_kw('world',	'STR', "universe")
	sect1.add_sect(subsect)

	sect2=getkw.Section('sect2')
	sect2.add_kw('foobar',	'INT', -42)
	sect2.add_kw('raboof',	'DBL', -74.0e-6)
	sect2.add_kw('oof',	'DATA')
	top.add_sect(sect2)
	
	if inpfil is None:
		inkw=getkw.Getkw(top) # we use ourselves as stencil
	else:
		print "Parsing '%s'" % inpfil
		input=getkw.GetkwParser()
		inkw=input.parseFile(inpfil)
		inkw.sanitize(top)

	topsect=inkw.get_topsect()

	if opts.debug:
		inkw.setkw('debug', opts.debug)
		dbg=opts.debug
	else:
		dbg = 1 # testing only

	inkw.run_callbacks(top)
	
	print "Running with the following input..."
	tmpfile='_example.inp'
	tfd=open(tmpfile,'w')
	print >> tfd, inkw.top
	tfd.close()
	os.system('cat < ' + tmpfile)
	if dbg == 0:
		os.unlink(tmpfile)

# Sample callback to check sanity. Using callbacks very sophisticated 
# inputs can be handled
def verify_top(top):
	warn="Warning: The '%s' option will be ignored in toplevel"
	err="Error: Required option '%s' not set for start!"
	required=()
	ignored=()

	for i in required:
		if not check_opt(top, i):
			print err % (i)
			sys.exit(1)
	for i in ignored:
		if check_opt(top, i):
			print warn % (i)
	
def check_opt(sect,key):
	try:
		k=sect[key]
	except:
		print 'You have a typo in the code for key', key
		sys.exit(1)
	if k is not None:
		if k.is_set():
			return True
	return False
		

if __name__ == '__main__':
	main()

