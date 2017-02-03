num_to_word_letters = 'ABCDEFGHIJKLMNOPQRTUVWXYZ'
num_to_word_history = {}
def num_to_word(num):
    if num not in num_to_word_history:
        base = len(num_to_word_letters)
        place = 1
        result = []
        while num >= 0:
            result.insert(0, num_to_word_letters[(num % base ** place) // base ** (place - 1)])
            num -= base ** place
            place += 1
        num_to_word_history[num] = ''.join(result)
    return num_to_word_history[num]

