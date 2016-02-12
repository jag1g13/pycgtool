class Section:
    __slots__ = ["name", "_lines"]

    def __init__(self, name=None):
        self.name = name
        self._lines = []

    def __len__(self):
        return len(self._lines)

    def __iter__(self):
        return iter(self._lines)

    def add_line(self, line):
        self._lines.append(line)


class DuplicateSectionError(Exception):
    pass


class CFG:
    __slots__ = ["filename", "_sections", "_section_names", "_iter_section"]

    def __init__(self, filename=None):
        self.filename = filename
        self._sections = {}
        self._section_names = []

        with open(self.filename) as f:
            curr_section = None
            for line in f:
                line = line.strip("\n")
                if line == "":
                    continue
                if line[0] in [";", "#"]:
                    continue
                if line[0] == "[":
                    curr_section = line.strip("[ ]")
                    if curr_section in self._sections:
                        raise DuplicateSectionError("Section {0} appears twice in file {1}".format(
                            curr_section, self.filename
                        ))
                    self._section_names.append(curr_section)
                    self._sections[curr_section] = Section(name=curr_section)
                    continue

                toks = line.split()
                self._sections[curr_section].add_line(toks)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def __len__(self):
        return len(self._section_names)

    def __iter__(self):
        self._iter_section = -1
        return self

    def __next__(self):
        self._iter_section += 1
        if self._iter_section >= len(self._section_names):
            raise StopIteration
        sec = self._section_names[self._iter_section]
        return self._sections[sec]

    def __contains__(self, item):
        return item in self._section_names

    def __getitem__(self, item):
        return self._sections[item]

