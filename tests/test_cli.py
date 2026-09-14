from protaine.cli import build_parser


def _parse_run_args(*extra):
    parser = build_parser()
    return parser.parse_args(["run", "--protein-fasta", "protein.fasta", *extra])


def test_structure_and_codon_weight_default_to_half_half():
    args = _parse_run_args()
    assert args.structure_weight == 0.5
    assert args.codon_weight == 0.5


def test_structure_and_codon_weight_are_configurable():
    args = _parse_run_args("--structure-weight", "0.8", "--codon-weight", "0.2")
    assert args.structure_weight == 0.8
    assert args.codon_weight == 0.2


def test_check_off_target_defaults_to_false():
    args = _parse_run_args()
    assert args.check_off_target is False


def test_check_off_target_flag_enables_it():
    args = _parse_run_args("--check-off-target")
    assert args.check_off_target is True


def test_predict_3d_defaults_to_disabled():
    args = _parse_run_args()
    assert args.predict_3d is False
    assert args.predict_3d_timeout == 300


def test_predict_3d_flag_and_timeout_are_configurable():
    args = _parse_run_args("--predict-3d", "--predict-3d-timeout", "60")
    assert args.predict_3d is True
    assert args.predict_3d_timeout == 60
