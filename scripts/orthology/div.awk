BEGIN { OFS="\t"}

{
    if ($4 != 0) {
        quotient = $3 / $4
        $15 = sprintf("%.3f", quotient)
    }
    else {
        $15 = "UND"
    }
    print
}
