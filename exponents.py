"""
Module exponents defines functions which enable translations between floats
and litteral expressions which are used to name files and directories
throughout this project.
These litteral expressions work as scientific notation and are read as follows
	e.g., k123456
	      ||`---'
	      ||  |_ decimals
	      ||____ integer part
	      |_____ exponent
	with a = -11 to z = 14, we then have k = -1, thus k123456 = 1.23456e-1,
	moreover, a capital letter preceding digits means the number is negative,
	e.g., K123455 = -1.23456e-1.

We set ONCE AND FOR ALL the number of significant figures with the variable
significant_figures.
"""

from numpy import floor, log10, abs, seterr
seterr(divide='ignore')	# ignore division by 0

# SIGNIFICANT FIGURES

significant_figures = 4	# number of significant figures displayed in file names

# TRANSLATIONS

_exponents = {'a':'-11', 'b':'-10', 'c':'-9', 'd':'-8', 'e':'-7', 'f':'-6',
	'g':'-5', 'h':'-4', 'i':'-3', 'j':'-2', 'k':'-1', 'l':'0', 'm':'1',
	'n':'2', 'o':'3', 'p':'4', 'q':'5', 'r':'6', 's':'7', 't':'8', 'u':'9',
	'v':'10', 'w':'11', 'x':'12', 'y':'13', 'z':'14'}	# associative array of exponents and corresponding letters

def letters_to_float(let):
	"""
	Translates litteral expression to float.

	Parameters
	----------
	let : string
		Litteral expression.

	Returns
	-------
	flo : float
		Corresponding float.
		NOTE: returns 0 if expression is incorrect.
	"""

	expr = str(let[1]						# integer part
		+ '.' + let[2:]						# decimals
		+ 'e' + _exponents[let[0].lower()]	# exponent
	)

	try:
		return (1 - 2*let[0].isupper())*float(expr)
	except ValueError: return 0

def float_to_letters(flo):
	"""
	Translates float to litteral expression.

	Parameters
	----------
	flo : float
		Corresponding float.

	Returns
	-------
	let : string
		Litteral expression.
		NOTE: returns l0000 (= 0) if expression is incorrect.
	"""

	try:
		expo_int = int(floor(log10(abs(flo))))			# integer exponent
		expo_let = list(_exponents.keys())[list(_exponents.values()).index(
			'%i' % expo_int								# corresponding letter exponent
			)]
	except (OverflowError, KeyError): return 'l0000'	# return 0 if exponent not attainable or zero
	if flo < 0: expo_let = expo_let.upper()				# make upper case if float is negative

	digi = int(abs(flo) * (10**(significant_figures - expo_int - 1)))	# digits in litteral expression

	return '%s%i' % (expo_let, digi)
