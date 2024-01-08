import sys
import itertools

def process_files(file_a, file_b):
    # read lines from file_b into a list
    with open(file_b, 'r') as f:
        b_lines = [line.strip().replace("\t","").replace(" ","") for line in f]

    # read lines from file_a, remove the last tab
    with open(file_a, 'r') as f:
        a_lines = [line.strip().replace("\t","").replace(" ","") for line in f]

    # check overlap and concat lines
    result = []
    for astr, bstr in itertools.combinations(a_lines, 2):  # get all pairs of lines
        for cstr in b_lines:
            # find the overlap between astr, cstr and bstr, cstr
            for i in range(1,len(cstr)):
                for j in range(i, len(cstr)):
                    if astr.endswith(cstr[:i]) and bstr.startswith(cstr[j:]):
                        concat_str = astr + cstr[i:j] + bstr[len(cstr[j:]):]
                        print(astr, cstr, bstr)
                        result.append(concat_str)
    return result

def main():
    # Check if correct arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python script.py <file_a> <file_b>")
        sys.exit(1)

    file_a = sys.argv[1]
    file_b = sys.argv[2]
    results = process_files(file_a, file_b)

    print("------------")
    for line in results:
        print(line)

if __name__ == "__main__":
    main()
