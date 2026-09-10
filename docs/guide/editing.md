# Moving and editing

Every change here is undoable: ++ctrl+z++ undoes, ++ctrl+y++ redoes, and the
history holds the last hundred edits.

## The manipulator

A manipulator sits on whatever is selected.

![sliding and turning a fragment with the manipulator](../assets/handles.gif)

| Drag | Moves the selection |
|---|---|
| a red, green or blue **arrow** | along x, y or z |
| a red, green or blue **ring** | round that axis |
| the **centre dot** | in the plane of the screen |

Hold ++shift++ while dragging to snap: 0.1 Å for a slide, 15° for a turn.
++g++ hides the handles when they are in the way.

Arrows win over rings where they overlap, and an axis seen almost end-on is not
offered at all — it would be too short on screen to drag accurately.

## Picking up a whole fragment

++f++ grows the selection to everything connected to it. Select one atom of a
ligand or a substituent, press ++f++, and the whole group is on the manipulator.

## An exact value

To set a distance, angle or dihedral to a number rather than drag it there,
select two, three or four atoms and press ++e++. A dialog adjusts the value with
a live preview, and everything beyond the last atom you selected moves with it
— the way GaussView does it — so select the atom on the side you want to move
last. Inside a ring there is no "beyond", and only that end atom moves.

## Deleting and rebuilding

- ++x++ deletes the selected atoms
- ++b++ recomputes the bonds, after a change large enough to make or break one

## Saving the result

++ctrl+s++ writes the edited structure as `xyz`, a Gaussian input (`.gjf`), an
ORCA input (`.inp`), a Q-Chem input (`.in`), `pdb`, or an MDL molfile (`.mol`) —
ready to submit.
