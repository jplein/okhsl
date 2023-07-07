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
	h := 283. / 360.
	s := 1.
	l := 10. / 100.

	actual := Okhsl_to_srgb(HSL{h: h, s: s, l: l})
	// Expected values come from the JS implementation which returns r, g, and b
	// between 0 and 255; our function, ported from C++, returns r, g, and b
	// between 0 and
	expected := []float64{21.025869038964053 / 255, 0.013077429697195551 / 255, 72.70843103040897 / 255}

	if !compare3Tuple([]float64{actual.r, actual.g, actual.b}, expected, 0.001) {
		t.Errorf("Unexpected value for 283 / 100 / 10: got '%s', expected '%s'", fmt.Sprint(actual), fmt.Sprint(expected))
	}
}

func TestSRGBToOKHSL(t *testing.T) {
	actual := Srgb_to_okhsl(RGB{21. / 255., 0, 73. / 255.})
	expected := []float64{0.78, 1.00, 0.10}

	if !compare3Tuple(
		[]float64{actual.h, actual.s, actual.l},
		expected,
		0.01,
	) {
		t.Errorf("Unexpected value for [141, 43, 16]: Expected %s but found %s", fmt.Sprint(expected), fmt.Sprint(actual))
	}
}
