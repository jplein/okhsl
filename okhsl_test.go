package okhsl

import (
	"fmt"
	"math"
	"testing"
)

func compareFloat(actual, expected, epsilon float64) bool {
	if actual == expected {
		return true
	}

	if math.Abs(actual-expected) < epsilon {
		return true
	}

	return false
}

func compare3Tuple(actual, expected []float64, epsilon float64) bool {
	aMatches := compareFloat(actual[0], expected[0], epsilon)
	bMatches := compareFloat(actual[1], expected[1], epsilon)
	cMatches := compareFloat(actual[2], expected[2], epsilon)

	return aMatches && bMatches && cMatches
}

func TestOKHSLToSRGB(t *testing.T) {
	actual := OKHSLToSRGB(283.0/360.0, 100.0/100.0, 10.0/100.0)
	expected := []float64{21.025869038964053, 0.013077429697195551, 72.70843103040897}

	if !compare3Tuple(actual, expected, 0.001) {
		t.Errorf("Unexpected value for 283 / 100 / 10: got '%s', expected '%s'", fmt.Sprint(actual), fmt.Sprint(expected))
	}
}

func TestSRGBToOKHSL(t *testing.T) {
	actual := SRGBToOKHSL(21, 0, 73)
	expected := []float64{0.78, 1.00, 0.10}

	if !compare3Tuple(actual, expected, 0.01) {
		t.Errorf("Unexpected value for [141, 43, 16]: Expected %s but found %s", fmt.Sprint(expected), fmt.Sprint(actual))
	}
}
