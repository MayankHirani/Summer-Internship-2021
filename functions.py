
import numpy as np
import math


# Create a one-hot matrix of a sequence
def onehot(sequence):

	one_hot_matrix = [ ]

	# print(sequence)

	for letter in sequence:
		# Each letter assigned to a different matrix
		dict1 = {'A':[1,0,0,0], 'C':[0,1,0,0], 'G':[0,0,1,0], 'T':[0,0,0,1], 'N':[0.25, 0.25, 0.25, 0.25]}

		# Add the matrix to the list for writing to the file
		one_hot_matrix.append(dict1[letter])

	return one_hot_matrix


# Turn a WTMX into a PWM numpy array
def pwm_from_wtmx(wtmx):

	# Create the PWM from each line
	PWM = [[float(x) for x in line.strip('\n').split('\t')] for line in wtmx[1:-1]]

	# Convert to numpy array
	PWM = np.array(PWM)

	return PWM


# Change the values to the result of the equation
def changeValues(pwm):

	# Change the values of the ones that are 0 to add white noise
	for row in range(len(pwm)):
		for i in range(len(pwm[row])):
			if pwm[row][i] == 0:
				pwm[row][i] = 10 ** -5

	# Change each value to the probability
	for row in range(len(pwm)):
		rowSum = sum(pwm[row])
		for i in range(len(pwm[row])):
			pwm[row][i] = pwm[row][i] / rowSum

	# Plug each probability into the equation
	for row in range(len(pwm)):
		for i in range(len(pwm[row])):
			# log2(val / bk), bk=0.25
			pwm[row][i] = math.log( pwm[row][i] / 0.25 , 2 )

	return pwm


# Calculate the score of a sequence given a PWM (Slower way of calculateBestScore())
def calculateScore(sequence, pwm):

	# Where each letter is located in the PWM
	index = {'A':0, 'C':1, 'G':2, 'T':3}

	score = 0

	# Go through the sequence and add up the scores at each index
	for i, letter in sequence:
		score += pwm[i][ index[letter] ]

	return score


# Get the maximum possible score of a PWM, to make sure that the values are less than the max
def calculateBestScoreOfPWM(pwm):

	score = 0

	for row in pwm:
		score += max(row)

	return abs(score)

# Calculate the best score of a sequence in onhot form, by sliding it along
def calculateBestScore(sequence, pwm):

	scores = []

	# Slide it along
	for i in range(0, len(sequence) - len(pwm) + 1):
		section = sequence[i : i + len(pwm)]

		multiplied = np.multiply(section, pwm)
		scores.append(np.sum(multiplied))

	bestScorePossible = calculateBestScoreOfPWM(pwm)
	maxScore = max(scores)

	# Check that the score isnt above the max possible
	if (maxScore - bestScorePossible) < (0.1 * bestScorePossible):
		# Return e to the power of the difference
		# Change this method to change the way it is scored
		return np.exp(maxScore - bestScorePossible)
		# return maxScore
	else:
		return "Error, score higher than max\nMax possible score: " + str(abs(calculateBestScoreOfPWM(pwm))) + "\nHighest Score: " + str(max(scores))







