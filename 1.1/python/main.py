def int_to_binary(integer):
    # 指数部分转化为二进制，左大右小，11位
    int_binary = []
    while integer >= 1:
        int_binary.append(int(integer % 2))
        integer = integer // 2
    while len(int_binary) < 11:
        int_binary.append(0)
    return list(reversed(int_binary))


def binary_to_int(binary_str):
    # 二进制指数转化为十进制
    integer, i = 0, 0
    while i < len(binary_str):
        integer += int(binary_str[-i-1])*2**i
        i += 1
    return integer


def fra_to_binary(fraction):
    # 小数部分转化为二进制，左大右小，53位
    fra_binary = []
    i = 0
    while i < 52:
        fraction *= 2
        if fraction >= 1:
            fraction -= 1
            fra_binary.append(1)
        else:
            fra_binary.append(0)
        i += 1
    return fra_binary


def binary_to_fra(binary_str):
    # 小数部分转化为十进制
    fra, i = 0, 0
    while i < len(binary_str):
        fra += int(binary_str[-i-1])*2**(-i-1)
        i += 1
    return fra


def move_decimal(number):
    # 对输入的十进制数做乘除2的操作使其在1-2之间，输出处理后的数字和移动的位数
    i = 0
    while number >= 2:
        number /= 2
        i += 1
    while number < 1:
        number *= 2
        i -= 1
    return number - 1, i


def decimal_to_binary_doublefloat(number):
    if number < 0:
        sign = 1
    else:
        sign = 0
    number = abs(number)
    mpart_d, epart_d = move_decimal(number)
    if epart_d >= 1024:
        print("Number overflow! Doublefloat number should be less than 1.79*10^308.")
        epart_b = int_to_binary(2047)
    elif epart_d <= -1023:
        print("Number too small!")
        epart_b = int_to_binary(0)
    else:
        epart_b = int_to_binary(epart_d + 1023)
    mpart_b = list(reversed(fra_to_binary(mpart_d)))
    return [sign], epart_b, mpart_b


def binary_doublefloat_to_decimal(number_str):
    if len(number_str) != 64:
        print("Not a legitimate binary double float: length is %s instead of 64" % len(number_str))
        raise Exception("Input error!")
    sign = int(number_str[0])
    epart_b = number_str[1:12]
    epart_d = binary_to_int(epart_b)-1023
    if epart_d == -1023:
        print("Number is considered infinite small!")
    elif epart_d == 1024:
        print("Number is considered infinite large!")
    mpart_b = number_str[12:64]
    mpart_d = binary_to_fra(mpart_b)+1
    return (-1)**sign*2**epart_d*mpart_d


if __name__ == '__main__':
    for i in decimal_to_binary_doublefloat(float(input("Input a decimal number: "))):
        for j in i:
            print(j, end='')
    print()
    print(binary_doublefloat_to_decimal(input("Input a binary double float number: ")))
