




def head():
    print('<head>')
    # print('<link rel="stylesheet" type="text/css" href="evn_calculator.css">')
    print('<style>\n.checkboxes label { white-space: nowrap; }</style>')
    print('</head>')


def selections(tag, title, keys, values):
    print('<div>')
    print('<label for="{0}">{1}</label><br>'.format(tag, title))
    print('<select name="{0}" id="{0}">'.format(tag))
    for key, value in zip(keys, values):
        print('<option value="{0}">{1}</option>'.format(key, value))

    print('</select><br><br>')
    print('</div>')


def checkboxes(tag, name, keys, values):
    print('<div class="{}">'.format(tag))
    for key, value in zip(keys, values):
        print('<label class="{0}"><input type="{0}" name="{1}" value="{2}">{3}</label>'.format(
              tag, name, key, value))

    print('</div>')



