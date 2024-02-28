//
// test-stack.cpp
//
// Copyright (c) 2020 Carlos Braga
// This program is free software; you can redistribute it and/or modify it
// under the terms of the MIT License. See accompanying LICENSE.md or
// https://opensource.org/licenses/MIT.
//

#include "catch2/catch.hpp"
#include "algods/linked/stack.h"

/// ---- Test Stack Data Structure ---------------------------------------------
///
long stack_calculate(std::string &S) {
    //
    // Its an error if we try to operate on an empty string
    //
    if (S.empty()) {
        return -1;
    }

    //
    // Construct a stack of operators and a stack of operands
    //
    Stack<long> stack_num;
    for (auto it = S.begin(); it != S.end(); ++it) {
        if (*it == '0' ||
            *it == '1' ||
            *it == '2' ||
            *it == '3' ||
            *it == '4' ||
            *it == '5' ||
            *it == '6' ||
            *it == '7' ||
            *it == '8' ||
            *it == '9') {
            stack_num.push(*it - '0');   // char 0
        } else if (*it == '+' || *it == '*') {
            if (stack_num.size() < 2) {
                return -1;
            }

            long long d1 = static_cast<long long>(stack_num.top());
            stack_num.pop();

            long long d2 = static_cast<long long>(stack_num.top());
            stack_num.pop();

            long long d3 = 0;
            if (*it == '+') {
                d3 = d1 + d2;
            } else if (*it == '*') {
                d3 = d1 * d2;
            }

            if (d3 > std::numeric_limits<long>::max()) {
                return -1;  // overflow error
            }

            stack_num.push(static_cast<long>(d3));
        } else {
            return -1;  // its not a digit nor an operator
        }
    }

    //
    // Digits stack is empty after we completed all operations.
    //
    if (stack_num.is_empty()) {
        return -1;
    }

    return stack_num.top();
}


TEST_CASE("Stack") {
    // Vector of strings
    std::vector<std::string> text {
        "Lorem",
        "ipsum",
        "dolor",
        "sit",
        "amet,",
        "consectetur",
        "adipiscing",
        "elit.",
        "Morbi",
        "sodales",
        "eleifend",
        "magna,",
        "in",
        "laoreet",
        "sapien."};

    SECTION("Test Stack API") {
        {
            //
            // Check initial state
            //
            Stack<std::string> stack_a;
            REQUIRE(stack_a.is_empty() == true);
            REQUIRE(stack_a.size() == 0);

            //
            // Check add()
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                stack_a.push(*it);
            }
            REQUIRE(stack_a.size() == text.size());

            //
            // Check stack iterator
            //
            auto it = stack_a.begin();
            auto jt = text.rbegin();
            while (it != stack_a.end() && jt != text.rend()) {
                REQUIRE(*it == *jt);
                ++it;
                ++jt;
            }

            //
            // Check stack clear()
            //
            stack_a.clear();
            REQUIRE(stack_a.is_empty() == true);
            REQUIRE(stack_a.size() == 0);
        }

        {
            //
            // Check initial state
            //
            Stack<std::string> stack_a;
            REQUIRE(stack_a.is_empty() == true);
            REQUIRE(stack_a.size() == 0);

            //
            // Check Stack copy constructors
            //
            for (auto it = text.begin(); it != text.end(); ++it) {
                stack_a.push(*it);
            }

            {
                Stack<std::string> stack_b(stack_a);
                REQUIRE(stack_b.size() == stack_a.size());
                auto it = stack_a.begin();
                auto jt = stack_b.begin();
                while (it != stack_a.end() && jt != stack_b.end()) {
                    REQUIRE(*it == *jt);
                    ++it;
                    ++jt;
                }
            }

            {
                Stack<std::string> stack_b = stack_a;
                REQUIRE(stack_b.size() == stack_a.size());
                auto it = stack_a.begin();
                auto jt = stack_b.begin();
                while (it != stack_a.end() && jt != stack_b.end()) {
                    REQUIRE(*it == *jt);
                    ++it;
                    ++jt;
                }
            }

            {
                Stack<std::string> stack_b;
                stack_b = stack_a;
                REQUIRE(stack_b.size() == stack_a.size());
                auto it = stack_a.begin();
                auto jt = stack_b.begin();
                while (it != stack_a.end() && jt != stack_b.end()) {
                    REQUIRE(*it == *jt);
                    ++it;
                    ++jt;
                }
            }
        }

        {
            //
            // Check Stack calculator
            //
            const long result = 21;
            std::string str_base("12+34+*");  // ((1+2) * (3+4)) = 21
            REQUIRE(stack_calculate(str_base) == 21);

            std::string str1(str_base);
            std::string str2;
            for (long i = 0; i < 8; ++i) {
                str1 += str_base + "+";
                str2 = str_base + std::to_string(i+2) + "*";
                REQUIRE(stack_calculate(str1) == (i+2)*result);
                REQUIRE(stack_calculate(str2) == (i+2)*result);
            }

            std::string str3("12*34*");
            const long count = 2048;
            for (long i = 0; i < count; ++i) {
                str3 += std::string("12*34*");
                std::string str4 = str3 + std::string((2*i+3),'+');
                REQUIRE(stack_calculate(str4) == (i+2)*(2+12));
            }
        }
    }
}

