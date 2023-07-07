package okhsl

import (
	"fmt"
	"math"
)

// h, s, l: hue, saturation, lightness values. h is in the range 0 to
// 360, and s and l in the range 0 to 100. An error is returned if any of the
// three is outside of these values.
//
// Returns a hex string with a leading #, e.g., "#ff00ff"
func OKHSLToSRGBHexString(h, s, l float64) (string, error) {
	if h < 0 || h > 360 {
		return "", fmt.Errorf("Expected h to be in the range 0 to 360 but found %f", h)
	}
	if s < 0 || s > 100 {
		return "", fmt.Errorf("Expected s to be in the range 0 to 100 but found %f", s)
	}
	if l < 0 || l > 100 {
		return "", fmt.Errorf("Expected l to be in the range 0 to 100 but found %f", l)
	}

	rgb := OKHSLToSRGB(float64(h)/360.0, float64(s)/100.0, float64(l)/100.0)

	r := int(math.Round(rgb[0]))
	g := int(math.Round(rgb[1]))
	b := int(math.Round(rgb[2]))

	if r < 0 {
		r = 0
	}
	if r > 255 {
		r = 255
	}

	if g < 0 {
		g = 0
	}
	if g > 255 {
		g = 255
	}

	if b < 0 {
		b = 0
	}
	if b > 255 {
		b = 255
	}

	c := r<<16 + g<<8 + b

	return fmt.Sprintf("#%06x", c), nil
}
