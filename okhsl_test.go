package okhsl_test

import (
	"fmt"
	"math"
	"testing"

	"github.com/jplein/okhsl"
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

func TestRGBHex(t *testing.T) {
	r := okhsl.RGB{R: 74., G: 99., B: 224.}
	hex, err := r.ToHex()

	if err != nil {
		t.Errorf("unexpected error converting to hex: %s", err.Error())
	}

	expected := "#4a63e0"

	if hex != expected {
		t.Errorf("unexpected hex value: Expected %s but got %s", expected, hex)
	}

	r = okhsl.RGB{R: 0, G: 0, B: 0}
	err = r.FromHex("#bb6a3e")

	if err != nil {
		t.Errorf("unexpected error reading from hex: %s", err.Error())
	}

	if r.R != 187. {
		t.Errorf("unexpected RGB value: expected R to be 187 but got %f", r.R)
	}

	if r.G != 106. {
		t.Errorf("unexpected RGB value: expected G to be 186 but got %f", r.G)
	}

	if r.B != 62. {
		t.Errorf("unexpected RGB value: expected B to be 62 but got %f", r.B)
	}
}

func TestHSLDataTypes(t *testing.T) {
	h := okhsl.HSLNormalized{240., 1., 0.5}
	err := h.Validate()

	if err == nil {
		t.Errorf("Expected an error but error is nil")
	}

	h = okhsl.HSLNormalized{0., 1., 0.5}
	err = h.Validate()

	if err != nil {
		t.Errorf("Expected no error but got %s", err.Error())
	}

	hsl := h.ToHSL()
	if hsl.H != 0. {
		t.Errorf("Invalid conversion of HSLNormalized to HSL: Expected %f for H but found %f", h.H/360, hsl.H)
	}

	if hsl.S != 100. {
		t.Errorf("Invalid conversion of HSLNormalized to HSL: Expected %f for S but found %f", h.S/100, hsl.S)
	}
}

func TestOKHSLToSRGB(t *testing.T) {
	h := 283. / 360.
	s := 1.
	l := 10. / 100.

	actual, err := okhsl.OKHSLToSRGBNormalized(okhsl.HSLNormalized{H: h, S: s, L: l})
	if err != nil {
		t.Errorf("unexpected error converting HSL to RGB: %s", err.Error())
	}

	// Expected values come from the JS implementation which returns r, g, and b
	// between 0 and 255; our function, ported from C++, returns r, g, and b
	// between 0 and
	expected := []float64{21.025869038964053 / 255, 0.013077429697195551 / 255, 72.70843103040897 / 255}

	if !compare3Tuple([]float64{actual.R, actual.G, actual.B}, expected, 0.001) {
		t.Errorf("Unexpected value for 283 / 100 / 10: got '%s', expected '%s'", fmt.Sprint(actual), fmt.Sprint(expected))
	}
}

func TestSRGBToOKHSL(t *testing.T) {
	actual, err := okhsl.SRGBToOKHSLNormalized(okhsl.RGBNormalized{R: 21. / 255., G: 0, B: 73. / 255.})
	if err != nil {
		t.Errorf("unexpected error converting RGB to HSL: %s", err.Error())
	}

	expected := []float64{0.78, 1.00, 0.10}

	if !compare3Tuple(
		[]float64{actual.H, actual.S, actual.L},
		expected,
		0.01,
	) {
		t.Errorf("Unexpected value for [141, 43, 16]: Expected %s but found %s", fmt.Sprint(expected), fmt.Sprint(actual))
	}
}
