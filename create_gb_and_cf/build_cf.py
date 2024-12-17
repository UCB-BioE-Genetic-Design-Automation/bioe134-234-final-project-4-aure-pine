from models.construction_file import *

def format_step(step):
    if isinstance(step, PCR):
        return f"PCR {step.forward_oligo} {step.reverse_oligo} {step.template} {step.output}"
    elif isinstance(step, GoldenGate):
        return f"GoldenGate {' '.join(step.dnas)} {step.enzyme} {step.output}"
    elif isinstance(step, Gibson):
        return f"Gibson {' '.join(step.dnas)} {step.output}"
    elif isinstance(step, Transform):
        return f"Transform {step.dna} {step.strain} {step.antibiotic} {step.output}"
    else:
        return f"# Unsupported step: {step.operation}"


###
# Different cloning strategies
###
def handle_golden_gate(components: dict, enzymes: list[str], plasmid: str):
    steps = [
        PCR(forward_oligo="FWD_primer", reverse_oligo="REV_primer", template=plasmid, output="PCR_product"),
        GoldenGate(dnas=["PCR_product"], enzyme=enzymes[0], output="GoldenGate_output")
    ]
    return steps


def handle_gibson(components: dict, enzymes: list[str], plasmid: str):
    steps = [
        PCR(forward_oligo="FWD_primer", reverse_oligo="REV_primer", template=plasmid, output="PCR_product"),
        Gibson(dnas=["PCR_product"], output="Gibson_output")
    ]
    return steps


def handle_restriction_ligation(components: dict, enzymes: list[str], plasmid: str):
    steps = [
        PCR(forward_oligo="FWD_primer", reverse_oligo="REV_primer", template=plasmid, output="PCR_product"),
        Digest(dna="PCR_product", enzymes=enzymes, fragSelect=0, output="Digested_output"),
        Ligate(dnas=["Digested_output", plasmid], output="Ligation_output")
    ]
    return steps

# Centralized dispatcher
STRATEGY_HANDLERS = {
    "GoldenGate": handle_golden_gate,
    "Gibson": handle_gibson,
    "RestrictionLigation": handle_restriction_ligation
}

def build_construction_file(components, cloning_strategy, enzymes, plasmid, antibiotic, output_filename: str):
    """Builds a Construction File describing experimental steps."""
    if "cds" not in components:
        raise ValueError("Components must include at least a 'cds' (coding sequence).")

    sequences = components.copy()

    # Dynamically fetch the strategy handler
    if cloning_strategy not in STRATEGY_HANDLERS:
        raise ValueError(f"Unsupported cloning strategy: {cloning_strategy}")

    steps = STRATEGY_HANDLERS[cloning_strategy](components, enzymes, plasmid)

    # Add the final transform step
    transform_output = steps[-1].output
    steps.append(Transform(dna=transform_output, strain="DH5alpha", antibiotic=antibiotic, output="Final_product"))

    # Write to the Construction File
    cf = ConstructionFile(steps=steps, sequences=sequences)
    cf_filename = f"{output_filename}.txt"
    with open(cf_filename, "w") as file:
        file.write("# Sequences\n")
        for name, seq in cf.sequences.items():
            file.write(f"{name}\t{seq}\n")
        file.write("\n# Steps\n")
        for step in cf.steps:
            file.write(format_step(step) + "\n")
    print(f"Construction File created: {cf_filename}")

    return cf_filename