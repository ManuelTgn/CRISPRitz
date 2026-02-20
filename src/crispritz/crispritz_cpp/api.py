from crispritz import _ternary_search_tree as tst  # type: ignore

def build_tree_cpp(sequence: str, pam: str, pam_length: int, pam_size: int, upstream: bool) -> None:
    tst.build_tree(sequence, pam, pam_length, pam_size, upstream)
