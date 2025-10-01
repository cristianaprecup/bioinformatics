#a DNA seq is fiven: S="ACGGGCATATGCGC".
# Make an app which is able to show the percentage of the components from the alphabet of the seq S.
# In otherwise the inpu of the seq S and the output is the alphabet of the seq and the percentage of each letter in the alphabet found in seq S

def calculate_percentage(sequence):
    total_length = len(sequence)

    letter_count = {}

    for letter in sequence:
        if letter in letter_count:
            letter_count[letter] += 1
        else:
            letter_count[letter] = 1

    letter_percentage = {
        letter: (count / total_length) * 100
        for letter, count in letter_count.items()
    }

    return letter_percentage


if __name__ == "__main__":
    S = "ACGGGCATATGCGC"

    percentages = calculate_percentage(S)
    print("Letter percentages:")
    for letter, percentage in percentages.items():
        print(f"{letter}: {percentage:.2f}%")