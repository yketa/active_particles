#! /bin/bash
#
# This bash shell script enables one to use functions
# active_particles.exponents.float_to_letters and
# active_particles.exponents.letters_to_float from the command line.
#
# Information about litteral notation of floats used in this project can be
# found in active_particles/exponents.py.

letters_to_float(){
	# Converts litteral expression to float expression.
	python -c "from active_particles.exponents import letters_to_float;\
	print(letters_to_float('$1'))"
}
export -f letters_to_float

float_to_letters(){
	# Converts float expression to litteral expression.
	python -c "from active_particles.exponents import float_to_letters;\
	print(float_to_letters($1))"
}
export -f float_to_letters
