# Style files

A style file is JSON. Save one from the [style editor](../guide/style.md), or
write one by hand — it only needs the keys you want to change, and the rest fall
back to the CYLview defaults below.

```json
{
  "name": "our group",
  "atom_scale": 0.42,
  "bond_color": "#1a1a1a",
  "palette": {"C": "#4d4d4d", "N": "#2f5fd0", "O": "#d13b2e"}
}
```

Use it from the viewer (the style editor's **Load…** button), the command line
(`--style ours.json`), or Python (`microscope.load_style("ours.json")`).

## Every key

This table is generated from the `Style` class itself, so it cannot fall out of
date.

<!-- style-reference -->

## Colours

A colour is `#rrggbb`, or three numbers from 0 to 1 such as `[0.3, 0.3, 0.3]`.
`null` switches off a colour that is optional — `"bond_color": null` splits each
bond by its two atoms instead of drawing it in one colour.

## The palette

`palette` maps element symbols to colours. Only the elements you list change;
anything else keeps its default, and an element nobody has coloured falls back
to its CPK colour.

```json
"palette": {"C": "#4d4d4d", "Pd": "#006985", "Ir": "#175487"}
```

## When a file is wrong

It names the problem instead of raising a traceback:

| File contains | Message |
|---|---|
| a key that is not part of a style | `not part of a style: colour` |
| an element that does not exist | `palette: 'Xx' is not an element` |
| a colour that is not `#rrggbb` | `background: 'red' is not #rrggbb` |
| a colour with the wrong count | `palette[C]: a colour needs three numbers` |
| broken JSON | `not valid JSON (…, line 1)` |
