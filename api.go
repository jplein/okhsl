package api

import "github.com/jplein/okhsl/types"

func PublicStub() types.RGB {
	r := 0.
	g := 255.
	b := 0.

	return types.RGB{
		R: r,
		G: g,
		B: b,
	}
}
