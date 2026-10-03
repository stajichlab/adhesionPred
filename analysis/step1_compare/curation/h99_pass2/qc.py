import html
import re


def norm(s):
    s = html.unescape(s)
    s = s.replace(" ", " ").replace("\xa0", " ")
    return re.sub(r"\s+", " ", s).strip()


def check(quote, fname):
    t = norm(open("texts/" + fname).read())
    q = norm(quote)
    return q in t
