# -*- coding: utf-8 -*-
from __future__ import annotations

import csv
import datetime as dt
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION_START
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[3]
REPORT_DIR = ROOT / "00_DELIVERABLES_ORGANIZED" / "01_CURRENT_LOG2_RESULTS" / "01_master_report"
OUT_DOCX = REPORT_DIR / "MASTER_GEMsplice_HumanGEM_iMAT_Isoform_Final_TR.docx"

FONT = "Calibri"
BLUE = "0F4A6E"
DARK = "172033"
MUTED = "5B6474"
LINE = "D7DEE8"
SOFT = "F7F9FC"
GREEN = "14532D"
AMBER = "8A4B0F"


def pth(*parts: str) -> Path:
    return ROOT.joinpath(*parts)


PATHS = {
    "old_master": REPORT_DIR / "MASTER_GEMsplice_HumanGEM_iMAT_Raporu_TR.pdf",
    "orig_plot1": pth("GEMsplice-master", "source code", "plot_outputs", "cancer_vs_normal", "pretty_png", "pretty_plot1.png"),
    "orig_plot2": pth("GEMsplice-master", "source code", "plot_outputs", "cancer_vs_normal", "pretty_png", "pretty_plot2.png"),
    "orig_plot3": pth("GEMsplice-master", "source code", "plot_outputs", "cancer_vs_normal", "pretty_png", "pretty_plot3.png"),
    "imat_plot1": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_hgemrules_tooltx2gene_ADobjectives", "plot_outputs", "plot1.png"),
    "imat_plot2": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_hgemrules_tooltx2gene_ADobjectives", "plot_outputs", "plot2.png"),
    "imat_plot3": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_hgemrules_tooltx2gene_ADobjectives", "plot_outputs", "plot3.png"),
    "relaxed_plot1": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_isoform_relaxed_mapping", "plot_outputs", "plot1.png"),
    "relaxed_plot2": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_isoform_relaxed_mapping", "plot_outputs", "plot2.png"),
    "relaxed_plot3": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_isoform_relaxed_mapping", "plot_outputs", "plot3.png"),
    "relaxed_summary": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_isoform_relaxed_mapping", "isoform_mapping_relaxed", "outputs", "mapping_summary.tsv"),
    "relaxed_group": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_isoform_relaxed_mapping", "outputs", "group_objective_summary.tsv"),
    "relaxed_samples": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_isoform_relaxed_mapping", "outputs", "sample_objectives.tsv"),
    "relaxed_objective": pth("imat_gemsplice_workspace", "gemsplice_imat_alz_control_isoform_relaxed_mapping", "isoform_mapping_relaxed", "outputs", "objective_reaction_mapping_check.tsv"),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fmt_num(value: str | float, nd: int = 2) -> str:
    try:
        x = float(value)
    except Exception:
        return str(value)
    if abs(x) >= 100:
        return f"{x:,.2f}".replace(",", "_").replace(".", ",").replace("_", ".")
    return f"{x:.2f}".replace(".", ",")


def set_run(run, size=10, bold=False, color=DARK, italic=False):
    run.font.name = FONT
    run._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = RGBColor.from_string(color)


def add_text(p, text, size=10, bold=False, color=DARK, italic=False):
    r = p.add_run(text)
    set_run(r, size=size, bold=bold, color=color, italic=italic)
    return r


def shade(cell, fill):
    pr = cell._tc.get_or_add_tcPr()
    shd = pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        pr.append(shd)
    shd.set(qn("w:fill"), fill)


def border(cell, color=LINE):
    pr = cell._tc.get_or_add_tcPr()
    borders = pr.first_child_found_in("w:tcBorders")
    if borders is None:
        borders = OxmlElement("w:tcBorders")
        pr.append(borders)
    for edge in ["top", "left", "bottom", "right"]:
        tag = "w:" + edge
        el = borders.find(qn(tag))
        if el is None:
            el = OxmlElement(tag)
            borders.append(el)
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), "4")
        el.set(qn("w:space"), "0")
        el.set(qn("w:color"), color)


def cell_text(cell, text, size=8.5, bold=False, fill=None, color=DARK):
    if fill:
        shade(cell, fill)
    border(cell)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(0)
    p.paragraph_format.line_spacing = 1.0
    add_text(p, str(text), size=size, bold=bold, color=color)


def section_title(doc, text, subtitle=None):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(3)
    add_text(p, text, size=18, bold=True, color=DARK)
    if subtitle:
        q = doc.add_paragraph()
        q.paragraph_format.space_after = Pt(8)
        add_text(q, subtitle, size=10, color=MUTED)


def callout(doc, label, text, fill="EEF6FF", border_color="BBD7F0", label_color=BLUE):
    table = doc.add_table(rows=1, cols=1)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.columns[0].width = Inches(7.15)
    cell = table.cell(0, 0)
    shade(cell, fill)
    border(cell, border_color)
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(3)
    p.paragraph_format.space_after = Pt(3)
    add_text(p, label + " ", size=9, bold=True, color=label_color)
    add_text(p, text, size=9, color=label_color)


def add_table(doc, headers, rows, widths=None, font_size=8):
    table = doc.add_table(rows=1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    if widths:
        for i, w in enumerate(widths):
            table.columns[i].width = Inches(w)
    for i, h in enumerate(headers):
        cell_text(table.cell(0, i), h, size=font_size, bold=True, fill="E8EDF3")
    for ridx, row in enumerate(rows, start=1):
        tr = table.add_row()
        fill = "F8FAFC" if ridx % 2 == 0 else None
        for cidx, val in enumerate(row):
            cell_text(tr.cells[cidx], val, size=font_size, fill=fill)
    return table


def code_block(doc, lines):
    table = doc.add_table(rows=1, cols=1)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    cell = table.cell(0, 0)
    shade(cell, "111827")
    border(cell, "111827")
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(4)
    for line in lines:
        add_text(p, line + "\n", size=8, color="F8FAFC")


def image(doc, path: Path, caption: str, width=6.9):
    if not path.exists():
        p = doc.add_paragraph()
        add_text(p, f"[Görsel bulunamadı: {path.name}]", size=9, color="8F1D1D")
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(path), width=Inches(width))
    c = doc.add_paragraph()
    c.alignment = WD_ALIGN_PARAGRAPH.CENTER
    c.paragraph_format.space_after = Pt(7)
    add_text(c, caption, size=8.2, color=MUTED, italic=True)


def page_break(doc):
    doc.add_page_break()


def setup_doc() -> Document:
    doc = Document()
    sec = doc.sections[0]
    sec.page_width = Inches(8.5)
    sec.page_height = Inches(11)
    sec.top_margin = Inches(0.55)
    sec.bottom_margin = Inches(0.55)
    sec.left_margin = Inches(0.62)
    sec.right_margin = Inches(0.62)

    styles = doc.styles
    styles["Normal"].font.name = FONT
    styles["Normal"]._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    styles["Normal"].font.size = Pt(10)
    styles["Normal"].font.color.rgb = RGBColor.from_string(DARK)
    for name, size in [("Heading 1", 18), ("Heading 2", 14), ("Heading 3", 11)]:
        st = styles[name]
        st.font.name = FONT
        st._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = RGBColor.from_string(DARK)
        st.paragraph_format.space_before = Pt(8)
        st.paragraph_format.space_after = Pt(3)
    return doc


def bullets(doc, items, size=9):
    for item in items:
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Inches(0.15)
        p.paragraph_format.space_after = Pt(1.5)
        add_text(p, "• ", size=size, color=BLUE)
        add_text(p, item, size=size, color=DARK)


def build_report():
    mapping_summary = {r["metric"]: r["value"] for r in read_tsv(PATHS["relaxed_summary"])}
    group_rows = read_tsv(PATHS["relaxed_group"])
    sample_rows = read_tsv(PATHS["relaxed_samples"])
    obj_rows = read_tsv(PATHS["relaxed_objective"])

    doc = setup_doc()

    # Cover
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    add_text(p, "MASTER TECHNICAL REPORT", size=11, bold=True, color=BLUE)
    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(6)
    add_text(p, "GemSplice → Human-GEM → iMAT → Isoform-Specific Mapping", size=26, bold=True)
    p = doc.add_paragraph()
    add_text(
        p,
        "Orijinal GemSplice mantığından başlayarak Human-GEM uyarlaması, tx2gene transcript mapping, log2(TPM+1), common iMAT modeli ve son relaxed isoform/reaction mapping denemesine kadar yapılan çalışmanın güncel master raporu.",
        size=11,
        color=MUTED,
    )
    callout(
        doc,
        "Kapsam:",
        "Bu rapor önceki master PDF okunarak güncellendi. Python denemelerinin genel validasyon amacı dışında final MATLAB tool mantığına dahil edilmeyen kısımları sadeleştirildi; en son relaxed isoform mapping toolu, çıktıları, plotları ve fallback raporu eklendi.",
    )
    add_table(
        doc,
        ["Bileşen", "Güncel Durum", "Ana Çıktı"],
        [
            ["Orijinal GemSplice", "Referans algoritma", "GPR min/max, patient loop, plot mantığı"],
            ["Human-GEM", "fbabreast yerine insan metabolik ağı", "Human-GEM GEMsplice giriş formatına uyarlandı"],
            ["Expression", "log2(TPM+1)", "sqrt kaldırıldı; TPM tek log ölçeğine indirildi"],
            ["iMAT", "3 AD + 3 Control common model", "12.931 reaction → 4.658 kept reaction"],
            ["Relaxed isoform mapping", "Yeni ayrı tool", "1.967 reaction evidence-supported, plotlar üretildi"],
        ],
        widths=[1.6, 2.6, 3.0],
        font_size=8.2,
    )
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    add_text(p, f"Generated: {dt.datetime.now().strftime('%d.%m.%Y %H:%M')}", size=8, color=MUTED)
    page_break(doc)

    # Executive summary
    section_title(doc, "1. Executive Overview", "Ne yaptık, ne değişti, ne korundu?")
    bullets(
        doc,
        [
            "Korunan çekirdek: GemSplice’in reaction expression değerlendirme mantığı, GPR içinde or=max / and=min yaklaşımı, hasta bazlı objective çözümü ve pathway plot yapısı.",
            "Değiştirilen ana katman: fbabreast/Recon yerine Human-GEM ve ardından expression temelli common iMAT-reduced Human-GEM kullanıldı.",
            "Expression ölçeği: TPM değerleri log2(TPM+1) ile dönüştürüldü; daha önce denenmiş sqrt(log2(TPM+1)) sıkıştırması kaldırıldı.",
            "Transcript mapping: ID_conversion yerine tx2gene/genes.tsv temelli ENST/ENSG bağlantısı kullanıldı.",
            "Son yeni çalışma: UniProt, Rhea ve KEGG annotation ile daha az konservatif relaxed isoform/reaction mapping üretildi ve ayrı yeni tool içinde GEMsplice çalıştırıldı.",
        ],
    )
    add_table(
        doc,
        ["Adım", "Amaç", "Sonuç"],
        [
            ["Original GemSplice", "Algoritma referansı", "3 plot ailesi ve GPR min/max mantığı referans alındı"],
            ["Human-GEM replacement", "İnsan metabolik modeli", "Human-GEM reaction/gene yapısı GEMsplice’e bağlandı"],
            ["tx2gene mapping", "ENST transcript verisini ENSG genlerine bağlama", "ENSG→ENST max yaklaşımıyla çalışır prototip"],
            ["iMAT", "Ağı biyolojik olarak kısıtlama", "Common 4.658 reaction model"],
            ["Relaxed isoform mapping", "Reaction-specific evidence artırma", "562 EC, 39 Rhea overlap, 1.366 gene-UniProt EC supported"],
        ],
        widths=[1.7, 2.6, 3.0],
        font_size=8.2,
    )
    page_break(doc)

    # Original GemSplice
    section_title(doc, "2. Orijinal GemSplice Referansı", "Model farklı olsa da korunması gereken algoritmik çekirdek.")
    bullets(
        doc,
        [
            "Orijinal GemSplice standart Recon’u çıplak haliyle kullanmıyor; fbabreast/Recon tabanlı model için yazarların hazırladığı transcript/isoform conversion tabloları var.",
            "Reaction expression içinde GPR mantığı korunuyor: or ilişkisi max(), and ilişkisi min() olarak değerlendiriliyor.",
            "Hasta verisinde eşleşmeyen gene/transcript için default 1 yaklaşımı var; no-gene reactionlar reaction_expression içinde boş kalabiliyor.",
        ],
    )
    add_table(
        doc,
        ["Orijinal Sayım", "Değer", "Yorum"],
        [
            ["Toplam reaction", "3.788", "fbabreast_trilevel modelinde"],
            ["Gene/transcript expression içeren reaction", "2.307", "reaction_expression token içeriyor"],
            ["Boş/no-gene reaction", "1.481", "expression ile bağlanmayan reactionlar"],
            ["Gene sayısı", "1.905", "tamamı transcript entry ile eşleşmiş"],
        ],
        widths=[2.4, 1.0, 4.0],
    )
    code_block(
        doc,
        [
            "Orijinal GPR mantığı:",
            "geneA or geneB  -> max(geneA, geneB)",
            "geneA and geneB -> min(geneA, geneB)",
            "",
            "evaluate_objective:",
            "reaction_expression -> numeric expression -> flux bounds -> trilevel FBA objectives",
        ],
    )
    image(doc, PATHS["orig_plot1"], "Orijinal GemSplice plot 1: objective değerlerinin grup bazlı dağılımı.", width=6.7)
    image(doc, PATHS["orig_plot2"], "Orijinal GemSplice plot 2: objective space ve Pareto/front çizgileri.", width=6.7)
    image(doc, PATHS["orig_plot3"], "Orijinal GemSplice plot 3: pathway-level relative flux change.", width=6.7)
    page_break(doc)

    section_title(doc, "3. Human-GEM Uyarlaması", "Aynı algoritmayı farklı model altyapısına taşıma.")
    bullets(
        doc,
        [
            "Human-GEM GPR kuralları çoğunlukla ENSG gene düzeyindedir; hasta Salmon verisi ise ENST transcript TPM düzeyindedir.",
            "İlk yaklaşımda her ENSG, tx2gene.tsv ile ENST transcript kümesine açıldı ve gene expression için max(ENST...) kullanıldı.",
            "Bu yaklaşım isoform-aware bir ara çözümdür; fakat hangi isoformun hangi reaction’ı taşıdığını reaction düzeyinde kanıtlamaz.",
        ],
    )
    code_block(
        doc,
        [
            "Mevcut ara çözüm:",
            "ENSG_X -> max(/ENST_A/, /ENST_B/, /ENST_C/)",
            "",
            "Hasta bazında:",
            "ENST_A=2.1, ENST_B=6.4, ENST_C=0.3",
            "gene expression = max(2.1, 6.4, 0.3) = 6.4",
        ],
    )
    callout(
        doc,
        "Düzeltme:",
        "Orijinal GemSplice isoformları ortalama almıyor; model/transcript mapping hazır olduğunda her model gene/isoform kendi expression girdisine bağlanıyor. Human-GEM’de bu bağlantıyı kurmak için ek annotation gerekir.",
        fill="FFF7ED",
        border_color="F0C48A",
        label_color=AMBER,
    )
    page_break(doc)

    # iMAT
    section_title(doc, "4. iMAT ile Common Human-GEM Modeli", "Düz Human-GEM yerine expression uyumlu reaction ağı.")
    add_table(
        doc,
        ["Karar", "Seçilen yaklaşım", "Neden"],
        [
            ["Model stratejisi", "Common iMAT model", "AD ve Control aynı reaction ağı üzerinde kıyaslansın"],
            ["Expression dönüşümü", "log2(TPM+1)", "TPM dağılımını tek ve okunabilir ölçeğe indirir"],
            ["High/low sınıflama", "25/75 persentil + top/bottom 100", "3+3 küçük dataset için konservatif seçim"],
            ["Solver", "intlinprog", "iMAT MILP problemi için mevcut MATLAB çözücü"],
            ["Protected reactions", "PK/LDH/GDH + biomass support", "GEMsplice objective reactionları düşmesin"],
        ],
        widths=[1.8, 2.4, 3.0],
    )
    code_block(
        doc,
        [
            "expr = log2(TPM + 1);",
            "lowTh  = prctile(expr(expr>0), 25);",
            "highTh = prctile(expr(expr>0), 75);",
            "tissueModel = iMAT_intlinprog_local(model, rxnExpr, lb, ub, coreRxns, protectedRxns);",
            "GEMsplice_RNAseq(tissueModel, patients);",
        ],
    )
    image(doc, PATHS["imat_plot1"], "iMAT/tx2gene AD-objective tool: sample objective space.", width=6.7)
    image(doc, PATHS["imat_plot2"], "iMAT/tx2gene AD-objective tool: objective pair plotları ve Pareto/front çizgileri.", width=6.7)
    page_break(doc)

    image(doc, PATHS["imat_plot3"], "iMAT/tx2gene AD-objective tool: pathway-level relative flux change.", width=6.8)
    callout(
        doc,
        "Yorum:",
        "Bu aşama final büyük veri analizi değildir; 3 AD + 3 Control ile yöntem/prototip doğrulamasıdır. Sonuçlar hipotez üretimi ve pipeline validasyonu düzeyinde yorumlanmalıdır.",
        fill="F0FDF4",
        border_color="B7E4C7",
        label_color=GREEN,
    )
    page_break(doc)

    # Isoform mapping plan and relaxed result
    section_title(doc, "5. Reaction-Specific Isoform Mapping", "Neden gerekli ve nasıl otomatikleştirildi?")
    bullets(
        doc,
        [
            "Human-GEM’de reaction genellikle ENSG düzeyinde yazılıdır; ENST/protein isoform hangi reaction’ı gerçekten taşıyor bilgisi doğrudan modelde yoktur.",
            "Bu yüzden genes.tsv içindeki ENSG→ENST→ENSP→UniProt zinciri, reactions.tsv içindeki Rhea/KEGG/MetaNetX reaction ID’leri ve UniProt/Rhea/KEGG EC annotationları birleştirildi.",
            "İlk strict pipeline yalnızca EC/Rhea güçlü eşleşmesini kabul etti; sonra relaxed pipeline KEGG EC, Rhea overlap ve gene-UniProt EC destekli mapping seviyelerini ayrı evidence olarak ekledi.",
        ],
    )
    add_table(
        doc,
        ["Kaynak dosya/API", "Verdiği bilgi", "Kullanım"],
        [
            ["genes.tsv", "ENSG, ENST, ENSP, UniProt", "Transcript/protein adayları"],
            ["reactions.tsv", "KEGG, Rhea, MetaNetX, Recon3D IDs", "Reaction dış annotation"],
            ["UniProt REST", "EC, catalytic activity, reviewed status", "Protein/enzyme kanıtı"],
            ["Rhea REST", "Reaction equation, EC, KEGG cross-ref", "Reaction EC kanıtı"],
            ["KEGG REST", "Reaction→enzyme EC", "Relaxed EC coverage artırımı"],
        ],
        widths=[1.7, 2.7, 3.0],
    )
    add_table(
        doc,
        ["Evidence/Fallback", "Reaction Sayısı", "Yorum"],
        [
            ["reaction_specific_ec_match", mapping_summary.get("reaction_specific_ec_match", "562"), "Güçlü reaction EC eşleşmesi"],
            ["reaction_specific_rhea_overlap", mapping_summary.get("reaction_specific_rhea_overlap", "39"), "UniProt catalytic Rhea ile reaction Rhea örtüşmesi"],
            ["gene_uniprot_ec_supported", mapping_summary.get("gene_uniprot_ec_supported", "1366"), "Gene/UniProt EC destekli orta seviye evidence"],
            ["protein_coding_enst_fallback", mapping_summary.get("protein_coding_enst_fallback", "994"), "Protein-coding transcript fallback"],
            ["default_1", mapping_summary.get("default_1", "1697"), "grRule/no-gene reaction; expression bağlanmadı"],
        ],
        widths=[2.4, 1.2, 4.0],
    )
    page_break(doc)

    section_title(doc, "6. Relaxed Isoform Mapping ile Çalıştırılan Yeni Tool", "Eski toollara dokunmadan yeni klasörde üretildi.")
    bullets(
        doc,
        [
            "Yeni tool klasörü: imat_gemsplice_workspace/gemsplice_imat_alz_control_isoform_relaxed_mapping",
            "GEMsplice akışı, evaluate_objective, flux solver ve plot mantığı aynı kaldı.",
            "Fark: reaction_expression artık relaxed isoform mapping candidate tablosundan geliyor.",
            "Fallback kullanılan reactionlar raporlandı; reaction-specific olmayanlar final biyolojik iddia olarak sunulmamalı.",
        ],
    )
    code_block(
        doc,
        [
            "Fallback sırası:",
            "1) Reaction-specific ENST / EC evidence varsa onu kullan",
            "2) Yoksa gene-UniProt EC supported mapping",
            "3) Yoksa protein-coding ENST fallback",
            "4) Yoksa GEMsplice default 1",
        ],
    )
    image(doc, PATHS["relaxed_plot1"], "Relaxed isoform mapping tool: 3 AD + 3 Control objective space.", width=6.7)
    page_break(doc)

    # Results
    section_title(doc, "7. Sonuçlar", "Relaxed isoform mapping tool objective ve pathway çıktıları.")
    group_table = []
    for r in group_rows:
        group_table.append(
            [
                r.get("group", ""),
                r.get("n", ""),
                fmt_num(r.get("mean_mitochondrial_complex_I", "")),
                fmt_num(r.get("mean_lactate_max", "")),
                fmt_num(r.get("mean_glutamate_lactate_max_min", "")),
            ]
        )
    add_table(
        doc,
        ["Grup", "n", "Mean Complex I", "Mean LDH max", "Mean GDH/LDH max-min"],
        group_table,
        widths=[1.3, 0.6, 1.6, 1.6, 2.0],
        font_size=8.3,
    )
    sample_table = []
    for r in sample_rows:
        sample_table.append(
            [
                r.get("sample", ""),
                r.get("group", ""),
                fmt_num(r.get("mitochondrial_complex_I", "")),
                fmt_num(r.get("lactate_max", "")),
                fmt_num(r.get("glutamate_lactate_max_min", "")),
            ]
        )
    add_table(
        doc,
        ["Sample", "Grup", "Complex I", "LDH max", "GDH/LDH max-min"],
        sample_table,
        widths=[1.6, 1.2, 1.4, 1.4, 1.9],
        font_size=7.8,
    )
    bullets(
        doc,
        [
            "Alzheimer grubunda mitochondrial Complex I ortalaması Control’den yüksek çıktı.",
            "LDH max Alzheimer’da hafif yüksek görünüyor.",
            "GDH/LDH max-min objective Alzheimer’da daha az negatif; bu küçük 3+3 dataset için ayrışma işareti olarak yorumlanabilir.",
        ],
        size=8.7,
    )
    page_break(doc)

    image(doc, PATHS["relaxed_plot2"], "Relaxed isoform mapping tool: objective pair plotları ve Pareto/front çizgileri.", width=6.8)
    page_break(doc)
    image(doc, PATHS["relaxed_plot3"], "Relaxed isoform mapping tool: pathway-level relative flux change.", width=6.8)
    bullets(
        doc,
        [
            "Pathway plotlarında lipid/fatty acid metabolism, arachidonic acid, omega-6 fatty acid metabolism ve heparan sulfate degradation gibi yollar dikkat çekiyor.",
            "Bu temalar Alzheimer literatürüyle biyolojik olarak tamamen kopuk değildir; ancak 3+3 küçük veri nedeniyle kesin sonuç değil, hipotez üretici bulgudur.",
            "Plot mantığı orijinal GemSplice ile aynıdır; sonuçların farklı çıkması model, iMAT ve yeni reaction_expression katmanından kaynaklanır.",
        ],
        size=8.7,
    )
    page_break(doc)

    # Objective mapping check
    section_title(doc, "8. Objective Reaction Mapping Kontrolü", "GEMsplice objective reactionları hangi evidence ile çalıştı?")
    obj_table = []
    for r in obj_rows:
        obj_table.append(
            [
                r.get("rxn", ""),
                r.get("objective", ""),
                r.get("mapping_source", ""),
                r.get("reaction_EC", ""),
                r.get("gene_UniProt_EC", "")[:60],
                "Evet" if r.get("needs_manual_review", "") == "1" else "Hayır",
            ]
        )
    add_table(
        doc,
        ["rxn", "Objective", "Mapping source", "Reaction EC", "Gene UniProt EC", "Manual review"],
        obj_table,
        widths=[1.0, 1.6, 1.8, 1.1, 1.4, 1.0],
        font_size=7.5,
    )
    callout(
        doc,
        "Önemli:",
        "LDH ve GDH objective reactionları reaction-specific EC evidence aldı. Mitochondrial Complex I için gene tarafında EC desteği var ama reaction-side EC/Rhea eşleşmesi yetersiz kaldığı için fallback/evidence-lower sınıfta değerlendirilmelidir.",
        fill="FFF7ED",
        border_color="F0C48A",
        label_color=AMBER,
    )

    section_title(doc, "9. Bilimsel Yorum ve Sınırlar", "Ne söyleyebiliriz, ne söylememeliyiz?")
    add_table(
        doc,
        ["Söylenebilir", "Söylenmemeli"],
        [
            ["Human-GEM/iMAT üzerinde GemSplice mantığı korunarak çalışan ayrı bir relaxed isoform mapping prototipi üretildi.", "Tüm reactionlar kesin isoform-specific çözüldü."],
            ["1.967 reaction EC/Rhea/UniProt evidence ile desteklendi.", "Fallback reactionlar reaction-specific kanıtlıdır."],
            ["3+3 dataset sonuçları yöntem doğrulama ve hipotez üretimi için değerlidir.", "Bu sonuçlar klinik/nihai Alzheimer sonucu olarak sunulabilir."],
            ["Plotlar orijinal GemSplice mantığıyla üretildi.", "Sonuçların orijinal meme kanseri çalışmasıyla aynı çıkması beklenir."],
        ],
        widths=[3.6, 3.6],
        font_size=8.2,
    )
    bullets(
        doc,
        [
            "En güçlü sonraki adım: daha büyük, doku-tutarlı AD/Control RNA-seq dataset ile iMAT modelini yeniden üretmek.",
            "Relaxed mapping pipeline aynı şekilde yeni büyük iMAT modeline uygulanabilir.",
            "Manual review_needed.tsv dosyası özellikle objective/pathway açısından önemli reactionlar için incelenmelidir.",
        ],
        size=8.8,
    )
    page_break(doc)

    section_title(doc, "10. Final Yol Haritası", "Bu çalışma nasıl tamamlanır?")
    add_table(
        doc,
        ["Aşama", "Ne yapılacak?", "Beklenen çıktı"],
        [
            ["A1", "Büyük AD/Control RNA-seq dataset seçimi", "Doku/batch bilgisi temiz sample metadata"],
            ["A2", "Salmon/TPM ve log2(TPM+1)", "ENST expression matrisi"],
            ["A3", "Common iMAT modeli", "Daha güvenilir curated Human-GEM"],
            ["A4", "Relaxed isoform mapping pipeline", "Evidence-tagged reaction_expression"],
            ["A5", "GEMsplice çalıştırma", "Objective/pathway plotlar"],
            ["A6", "Manual review", "Kritik fallback reactionların düzeltilmesi"],
        ],
        widths=[0.8, 3.6, 3.0],
    )
    callout(
        doc,
        "Son cümle:",
        "Bu rapordaki son tool, önceki Human-GEM uyarlamasından daha güçlüdür: GemSplice algoritmasını korur, iMAT ile ağı kısıtlar, UniProt/Rhea/KEGG evidence ile reaction_expression katmanını geliştirir ve fallbackleri açıkça raporlar.",
        fill="F0FDF4",
        border_color="B7E4C7",
        label_color=GREEN,
    )

    # Footer
    for section in doc.sections:
        footer = section.footer.paragraphs[0]
        footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
        add_text(footer, "Burak Keskin | GemSplice - Human-GEM - iMAT - Isoform Mapping Master Report", size=7, color=MUTED)

    doc.save(OUT_DOCX)
    return OUT_DOCX


if __name__ == "__main__":
    path = build_report()
    print(path)
