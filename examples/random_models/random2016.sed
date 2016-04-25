#! /usr/bin/sed -f
/^Warning: Label specification/ {
				g
				s/^@//
				p
				# If the hold buffer was empty, set it to ``@''
				# to mark that titles should be printed.
				s/^$/@/
				# but if it contains at least two characters, empty it.
				s/...*//
				h
				d
			}
$			g

/^-- benchmark: /	d
/^States  *[0-9]* *$/	d
/^Trans.  *[0-9]* *$/	d
/^Actions  *[0-9]* *$/	d
/^ *killme.pa[:. ]*$/	d
/^<49>/			d
/All states will have the same label in this model/ d
/^$/			d
/^----------------/	d
# lines that end with numbers, possibly with a dot: these lines are to be joined together (and printed at the end)
/ [0-9][.0-9]* *$/	{
				# trim trailing spaces:
				s/ *$//
				# prepend the hold buffer:
				H
				g
				# If the buffer starts with an @, then we append the word
				# to the first line of the buffer; otherwise, we delete it.
				# The number is added to the second line of the buffer.
				# If the first line is not empty and starts with @:
				s/^@\(..*\)\n\(.*\)\n\(.*[^ ]\)   *\([0-9.]*\)/@\1,\3\
\2,\4/
				# If the first line is empty except for an @:
				s/^@\n\(.*[^ ]\)   *\([0-9.]*\)/@\1\
\2/
				# If the buffer does not start with @ and is not empty:
				s/^\([^@].*\)\n.*[^ ]   *\([0-9.]*\)/\1,\2/
				# If the first line is completely empty:
				s/\n.*[^ ]   *\([0-9.]*\)/\1/
				h
				d
			}
