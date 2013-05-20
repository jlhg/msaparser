# msaparser - Parsing clustal file for sequence variation analysis
#
# Copyright (C) 2013, Jian-Long Huang
# Licensed under The MIT License
# http://opensource.org/licenses/MIT
#
# Author: Jian-Long Huang (jianlong@ntu.edu.tw)
# Version: 1.1
# Created: 2013.5.18


class Parser(object):

    def __init__(self, blockstar=0.9, blocklen=10, checknum=5):
        self._blockstar = float(blockstar)
        self._blocklen = int(blocklen)
        self._checknum = int(checknum)
        self._parse = 0

    def parse(self, handle, seqtype):
        """
        handle: clustal file
        alignment order: ss, rs, rc
        seqtype: 'n' for nucleotide; 'a' for amino acid
        """
        self._parse = 1
        if seqtype not in ('n', 'a'):
            raise Exception("Argument 'seqtype' must be 'n' or 'a'.")

        self._unknown_base = {
            'n': 'N',
            'a': 'X',
        }
        self.htmltag = {
            'clustal': 'alignment',
            'clutitle': 'title',
            'ss_name': 'ss',
            'rs_name': 'rs',
            'rc_name': 'rc',
            'ast': 'ast',
            'gap': 'gap',
            'dot': 'dot',
            'col': 'col',
            'count': 'count',
            'origin': 'o',
            'mutation': 'm',
            'seqtype': seqtype,
        }
        self._title = _Clutitle(tag=self.htmltag.get('clutitle'))
        self._ss = _CluSequence(tag_name=self.htmltag.get('ss_name'),
                                tag_gap=self.htmltag.get('gap'),
                                tag_seqtype=self.htmltag.get('seqtype'),
                                tag_var=self.htmltag.get('origin'))
        self._rs = _CluSequence(tag_name=self.htmltag.get('rs_name'),
                                tag_gap=self.htmltag.get('gap'),
                                tag_seqtype=self.htmltag.get('seqtype'),
                                tag_var=self.htmltag.get('mutation'))
        self._rc = _CluSequence(tag_name=self.htmltag.get('rc_name'),
                                tag_gap=self.htmltag.get('gap'),
                                tag_seqtype=self.htmltag.get('seqtype'),
                                tag_var=self.htmltag.get('origin'))
        self._star = _CluStar(tag_ast=self.htmltag.get('ast'),
                              tag_dot=self.htmltag.get('dot'),
                              tag_col=self.htmltag.get('col'),
                              tag_count=self.htmltag.get('count'))
        self.mutnum = 0
        self.mutposlist = []
        self.mutprofile = []
        self.block = []

        with open(handle, 'r') as fi:
            title = fi.readline().rstrip('\n')
            if title[0:7] == 'CLUSTAL':
                self._title.parse(title)
            else:
                raise Exception('Need a clustal file.')

            while True:
                line = fi.readline()
                if line == '':
                    # end of file
                    break
                elif line == '\n':
                    continue
                else:
                    self._ss.parse(line)
                    self._rs.parse(fi.readline())
                    self._rc.parse(fi.readline())
                    self._star.parse(fi.readline(), self._ss.linelength[-1])

        nogap_indexs = []
        for i in range(self._ss.get_sequence_length()):
            if self._ss._base[i] != '-' and self._rs._base[i] != '-' and self._rc._base[i] != '-':
                nogap_indexs.append(i)

        block = group_continuous_number(nogap_indexs)

        for i in block:
            real_block_length = i[1] - i[0] + 1
            real_starpct = self._star.get_starnum(i[0], i[1]) / float(real_block_length)

            if self._blocklen > real_block_length or self._blockstar > real_starpct:
                # does not pass
                continue

            self.block.append(str(i[0] + 1) + '..' + str(i[1] + 1) +
                              ' L=' + str(i[1] - i[0] + 1) +
                              ' SP=' + str(round(real_starpct, 2)))

            for j in range(i[0], i[1] + 1):
                base_ss = self._ss.get_base(j)
                base_rs = self._rs.get_base(j)
                base_rc = self._rc.get_base(j)

                if any([base_ss, base_rs, base_rc]) == self._unknown_base.get(seqtype):
                    continue

                if base_ss == base_rc and base_ss != base_rs:
                    # Sequence variation
                    if self._star.neighbor_star_check(j, self._checknum):
                        self.mutnum += 1
                        self.mutprofile.append('%s:%s' % (self._ss.get_bases(j - 2, j + 2),
                                                          self._rs.get_bases(j - 2, j + 2)))
                        self.mutposlist.append('%s-%s-%s' % (self._ss.get_raw_position(j),
                                                             self._rs.get_raw_position(j),
                                                             self._rc.get_raw_position(j)))
                        self._ss.add_mutposition(j)
                        self._rs.add_mutposition(j)
                        self._rc.add_mutposition(j)

    def get_clustal_html(self):
        if not self._parse:
            return None

        html = []
        html.append('<div class="%s">' % (self.htmltag.get('clustal')))
        html.append(self._title.get_html())
        html.append('<br><br>')

        for i in range(1, len(self._star.count)):
            html.append(self._ss.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append(self._rs.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append(self._rc.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append(self._star.get_html(self._star.count[i - 1], self._star.count[i]))
            html.append('<br>')

        html.append('</div>')

        return ''.join(html)


class _Clutitle(object):

    def __init__(self, tag='alignment'):
        self.line = None
        self.tag = tag

    def parse(self, line):
        self.line = line

    def get_html(self):
        return '<span class="%s">%s</span></br>' % (self.tag, self.line)


class _CluSequence(object):

    def __init__(self, tag_name, tag_gap, tag_seqtype, tag_var):
        import re
        self.name = None
        self.blank = None
        self.linelength = []
        self._base = []
        self.pattern = re.compile('(\S+)(\s+)(\S+)')
        self.tag_name = tag_name
        self.tag_gap = tag_gap
        self.tag_seqtype = tag_seqtype
        self.tag_var = tag_var
        self._mutposlist = []

    def parse(self, line):
        match = self.pattern.match(line)

        if match is not None:
            self.name = match.group(1)
            self.blank = match.group(2)
            self.linelength.append(len(match.group(3)))
            self._base = self._base + list(match.group(3))
        else:
            raise Exception('Format error.')

    def get_sequence_length(self):
        return len(self._base)

    def get_base(self, position):
        return self._base[position]

    def get_bases(self, start_position, stop_position):
        return ''.join(self._base[start_position:stop_position + 1])

    def get_raw_position(self, position):
        return len(''.join(self._base[0:position + 1]).replace('-', ''))

    def add_mutposition(self, position):
        self._mutposlist.append(position)

    def get_html(self, start_position, end_position):
        values = []
        values.append('<span class="%s">%s</span>' % (self.tag_name, self.name))
        values.append('&nbsp;' * len(self.blank))

        for i in range(start_position, end_position):
            if i in self._mutposlist:
                values.append('<span class="%s %s">%s</span>' % (self.tag_var, self._base[i], self._base[i]))
            elif i == '-':
                values.append('<span class="%s">%s</span>' % (self.tag_gap, self._base[i]))
            else:
                values.append('<span class="%s %s">%s</span>' % (self.tag_seqtype, self._base[i], self._base[i]))

        values.append('<br>')

        return ''.join(values)


class _CluStar(object):

    def __init__(self, tag_ast, tag_dot, tag_col, tag_count):
        self.blank = None
        self._base = []
        self.tag_ast = tag_ast
        self.tag_dot = tag_dot
        self.tag_col = tag_col
        self.tag_count = tag_count
        self.count = [0]

    def parse(self, line, linelength):
        line = line.rstrip('\n')
        self.blank = line[0: len(line) - linelength]
        self._base = self._base + list(line[len(line) - linelength:len(line)])

        if self.count:
            self.count.append(self.count[-1] + linelength)
        else:
            self.count.append(linelength)

    def get_starnum(self, start, end):
        star_num = 0

        for i in range(start, end + 1):
            if self._base[i] == '*':
                star_num += 1

        return star_num

    def neighbor_star_check(self, position, star_check_num):
        if position < star_check_num or position + star_check_num > len(self._base) - 1:
            # Mutation position is on sides, return False.
            return False

        # Check two sides
        star_num = 0

        for i in range(position + 1, position + 1 + star_check_num):
            if self._base[i] == '*':
                star_num += 1

        if star_num < star_check_num * 0.8:
            return False

        star_num = 0

        for i in range(position - star_check_num, position):
            if self._base[i] == '*':
                star_num += 1

        if star_num >= star_check_num * 0.8:
            return True
        else:
            return False

    def get_html(self, start_position, end_position):
        values = []
        values.append('&nbsp;' * len(self.blank))

        for i in range(start_position, end_position):
            if self._base[i] == '*':
                values.append('<span class="%s">%s</span>' % (self.tag_ast, self._base[i]))
            elif self._base[i] == ':':
                values.append('<span class="%s">%s</span>' % (self.tag_col, self._base[i]))
            elif self._base[i] == '.':
                values.append('<span class="%s">%s</span>' % (self.tag_dot, self._base[i]))
            elif self._base[i] == ' ':
                values.append('&nbsp;')
            else:
                raise Exception('Extra symbol was found.')

        values.append('<span class="%s">%s</span>' % (self.tag_count, end_position))
        values.append('<br>')

        return ''.join(values)


def group_continuous_number(numbers):
    if not numbers:
        yield 0, 0
    else:
        # yield a numbers of tuple of continuous numbers
        first = last = numbers[0]
        for n in numbers[1:]:
            if n - 1 == last:
                # part of the group
                last = n
            else:
                # not part of the group
                yield first, last
                first = last = n
        # yield the last group
        yield first, last
