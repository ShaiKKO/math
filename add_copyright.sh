#!/bin/bash

# Copyright header to add
COPYRIGHT="// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

"

# Function to add copyright if not present
add_copyright() {
    local file="$1"
    # Check if file already has copyright
    if ! grep -q "Copyright 2025 Colin MacRitchie" "$file"; then
        # Create temp file with copyright + original content
        echo "$COPYRIGHT" > "$file.tmp"
        cat "$file" >> "$file.tmp"
        mv "$file.tmp" "$file"
        echo "Added copyright to $file"
    else
        echo "Copyright already present in $file"
    fi
}

# Process all cubature headers
for file in include/boost/math/cubature/*.hpp include/boost/math/cubature/detail/*.hpp; do
    if [ -f "$file" ]; then
        add_copyright "$file"
    fi
done

echo "Copyright headers added to all cubature files"