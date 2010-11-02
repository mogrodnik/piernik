#!/usr/bin/python

import qa2,re,numpy

have_use           = re.compile("^\s{1,12}use\s")
remove_warn        = re.compile('''(?!.*QA_WARN .+)''', re.VERBOSE)
unwanted           = re.compile("(\s|&|\n)", re.VERBOSE)

def do_magic(files,options):
   name = files[0]
   glob  = []
   temp  = []
   for f in files[1:]:
      lines = open(f,'r').readlines()
      temp  = qa2.remove_amp(filter(remove_warn.match,lines),True)
      uses  = [f for f in filter(have_use.search, temp) if (re.match("\s{0,9}use "+name,f))]
      for f in uses:
         glob.extend(f.split("only: ")[1].strip().split(','))
   my_priv = numpy.unique([unwanted.sub('',f) for f in glob])
   print my_priv


if __name__ == "__main__":
   from optparse import OptionParser
   usage = "usage: %prog module_name FILES"
   parser = OptionParser(usage=usage)
   parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="make lots of noise [default]")
   parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose",
                  help="be vewwy quiet (I'm hunting wabbits)")
   parser.add_option("-f", "--force",
                  action="store_true", dest="force",
                  help="commit despite errors (It will be logged)")
   (options, args) = parser.parse_args()
   if len(args) < 1:
      parser.error("incorrect number of arguments")
   do_magic(args,options)

